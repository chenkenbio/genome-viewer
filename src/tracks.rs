use std::collections::HashMap;

use anyhow::{Context, Result, bail};
use bigtools::{BedEntry, BigBedRead, BigWigRead};
use http_range_client::HttpReader;

use crate::config::{
    AppConfig, TrackConfig, TrackKind, expand_path, is_remote_path, read_source_to_string,
};
use crate::model::{Feature, SignalBin, WindowFunction};

#[derive(Clone, Debug)]
pub struct TextTrack {
    pub features_by_chrom: HashMap<String, Vec<Feature>>,
}

#[derive(Clone, Copy, Default)]
struct AggregateBin {
    weighted_sum: f64,
    covered_bases: u64,
    max_value: f64,
    min_value: f64,
    observed_count: u64,
}

pub fn build_text_tracks(config: &AppConfig) -> Result<HashMap<String, TextTrack>> {
    let mut tracks = HashMap::new();

    for track in &config.tracks {
        if matches!(track.kind, TrackKind::Bed | TrackKind::Gtf) {
            let parsed = load_text_track(track)
                .with_context(|| format!("failed to preload text track '{}'", track.id))?;
            tracks.insert(track.id.clone(), parsed);
        }
    }

    Ok(tracks)
}

pub fn load_text_track(track: &TrackConfig) -> Result<TextTrack> {
    let text = read_source_to_string(&track.source)?;
    let mut features_by_chrom: HashMap<String, Vec<Feature>> = HashMap::new();

    for (line_number, line) in text.lines().enumerate() {
        let parsed = match track.kind {
            TrackKind::Bed => parse_bed_line(line),
            TrackKind::Gtf => parse_gtf_line(line),
            _ => None,
        };

        if let Some(feature) = parsed {
            features_by_chrom
                .entry(feature.chrom.clone())
                .or_default()
                .push(feature);
        } else if !skip_line(line) {
            tracing::debug!(
                track_id = %track.id,
                line_number = line_number + 1,
                "skipping unparseable line"
            );
        }
    }

    for features in features_by_chrom.values_mut() {
        features.sort_by_key(|feature| (feature.start, feature.end));
    }

    Ok(TextTrack { features_by_chrom })
}

pub fn infer_track_kind(source: &str) -> Option<TrackKind> {
    let lower = source.to_ascii_lowercase();
    if lower.ends_with(".bigwig") || lower.ends_with(".bw") {
        Some(TrackKind::BigWig)
    } else if lower.ends_with(".bigbed") || lower.ends_with(".bb") {
        Some(TrackKind::BigBed)
    } else if lower.ends_with(".gtf") || lower.ends_with(".gtf.gz") {
        Some(TrackKind::Gtf)
    } else if lower.ends_with(".bed") || lower.ends_with(".bed.gz") {
        Some(TrackKind::Bed)
    } else {
        None
    }
}

/// Returns true for all file extensions that igv.js can handle, including index files.
/// Used by the file browser to decide which files to show.
pub fn is_genomic_file(source: &str) -> bool {
    let lower = source.to_ascii_lowercase();
    // Data files
    lower.ends_with(".bam")
        || lower.ends_with(".cram")
        || lower.ends_with(".vcf.gz")
        || lower.ends_with(".vcf")
        || lower.ends_with(".bw")
        || lower.ends_with(".bigwig")
        || lower.ends_with(".bb")
        || lower.ends_with(".bigbed")
        || lower.ends_with(".bed")
        || lower.ends_with(".bed.gz")
        || lower.ends_with(".gtf")
        || lower.ends_with(".gtf.gz")
        || lower.ends_with(".gff")
        || lower.ends_with(".gff.gz")
        || lower.ends_with(".gff3")
        || lower.ends_with(".gff3.gz")
        || lower.ends_with(".seg")
        || lower.ends_with(".seg.gz")
        || lower.ends_with(".wig")
        || lower.ends_with(".bedgraph")
        || lower.ends_with(".bedgraph.gz")
        // Index files
        || lower.ends_with(".bai")
        || lower.ends_with(".crai")
        || lower.ends_with(".tbi")
        || lower.ends_with(".csi")
        || lower.ends_with(".idx")
}

/// Maps file extension to the igv.js format name.
/// Returns None for unrecognized extensions.
pub fn igvjs_format(source: &str) -> Option<&'static str> {
    let lower = source.to_ascii_lowercase();
    if lower.ends_with(".bam") {
        Some("bam")
    } else if lower.ends_with(".cram") {
        Some("cram")
    } else if lower.ends_with(".vcf.gz") || lower.ends_with(".vcf") {
        Some("vcf")
    } else if lower.ends_with(".bw") || lower.ends_with(".bigwig") {
        Some("bigwig")
    } else if lower.ends_with(".bb") || lower.ends_with(".bigbed") {
        Some("bigbed")
    } else if lower.ends_with(".bed") || lower.ends_with(".bed.gz") {
        Some("bed")
    } else if lower.ends_with(".gtf") || lower.ends_with(".gtf.gz") {
        Some("gtf")
    } else if lower.ends_with(".gff") || lower.ends_with(".gff.gz") {
        Some("gff")
    } else if lower.ends_with(".gff3") || lower.ends_with(".gff3.gz") {
        Some("gff3")
    } else if lower.ends_with(".seg") || lower.ends_with(".seg.gz") {
        Some("seg")
    } else if lower.ends_with(".wig") {
        Some("wig")
    } else if lower.ends_with(".bedgraph") || lower.ends_with(".bedgraph.gz") {
        Some("bedgraph")
    } else {
        None
    }
}

pub fn query_text_track(
    track: &TextTrack,
    chrom: &str,
    start: u64,
    end: u64,
    limit: usize,
) -> Vec<Feature> {
    let Some(features) = track.features_by_chrom.get(chrom) else {
        return Vec::new();
    };

    let start_index = features.partition_point(|feature| feature.end <= start);
    let mut matches = Vec::new();

    for feature in &features[start_index..] {
        if feature.start >= end {
            break;
        }
        if feature.end > start {
            matches.push(feature.clone());
        }
        if matches.len() >= limit {
            break;
        }
    }

    matches
}

pub fn query_bigwig(
    source: &str,
    chrom: &str,
    start: u64,
    end: u64,
    bins: usize,
    window_function: WindowFunction,
) -> Result<Vec<SignalBin>> {
    let start = u32::try_from(start).context("start exceeds BigWig coordinate range")?;
    let end = u32::try_from(end).context("end exceeds BigWig coordinate range")?;
    if end <= start {
        bail!("end must be greater than start");
    }

    let bins = bins.clamp(1, 4000);
    let width = u64::from(end - start);
    let bin_width = width.div_ceil(bins as u64).max(1);

    let mut aggregates = vec![
        AggregateBin {
            min_value: f64::INFINITY,
            ..AggregateBin::default()
        };
        bins
    ];

    if is_remote_path(source) {
        let mut reader = BigWigRead::open(HttpReader::new(source))
            .with_context(|| format!("failed to open remote BigWig '{}'", source))?
            .cached();
        let iter = reader.get_interval(chrom, start, end).with_context(|| {
            format!(
                "failed to query BigWig interval {}:{}-{}",
                chrom, start, end
            )
        })?;
        accumulate_bigwig_records(iter, &mut aggregates, start, end, bin_width, bins)?;
    } else {
        let path = expand_path(source)?;
        let mut reader = BigWigRead::open_file(&path)
            .with_context(|| format!("failed to open BigWig '{}'", path.display()))?
            .cached();
        let iter = reader.get_interval(chrom, start, end).with_context(|| {
            format!(
                "failed to query BigWig interval {}:{}-{}",
                chrom, start, end
            )
        })?;
        accumulate_bigwig_records(iter, &mut aggregates, start, end, bin_width, bins)?;
    }

    let response = aggregates
        .into_iter()
        .enumerate()
        .map(|(index, aggregate)| {
            let bin_start = u64::from(start) + (index as u64 * bin_width);
            let bin_end = (bin_start + bin_width).min(u64::from(end));
            let value = match window_function {
                WindowFunction::None | WindowFunction::Mean => {
                    if aggregate.covered_bases > 0 {
                        aggregate.weighted_sum / aggregate.covered_bases as f64
                    } else {
                        0.0
                    }
                }
                WindowFunction::Max => {
                    if aggregate.observed_count > 0 {
                        aggregate.max_value
                    } else {
                        0.0
                    }
                }
                WindowFunction::Min => {
                    if aggregate.observed_count > 0 {
                        aggregate.min_value
                    } else {
                        0.0
                    }
                }
                WindowFunction::Count => aggregate.observed_count as f64,
                WindowFunction::Density => {
                    if bin_end > bin_start {
                        aggregate.covered_bases as f64 / (bin_end - bin_start) as f64
                    } else {
                        0.0
                    }
                }
            };
            SignalBin {
                start: bin_start,
                end: bin_end,
                value,
            }
        })
        .collect();

    Ok(response)
}

pub fn query_bigbed(
    source: &str,
    chrom: &str,
    start: u64,
    end: u64,
    limit: usize,
) -> Result<Vec<Feature>> {
    let start = u32::try_from(start).context("start exceeds BigBed coordinate range")?;
    let end = u32::try_from(end).context("end exceeds BigBed coordinate range")?;
    if end <= start {
        bail!("end must be greater than start");
    }

    let mut features = Vec::new();
    if is_remote_path(source) {
        let mut reader = BigBedRead::open(HttpReader::new(source))
            .with_context(|| format!("failed to open remote BigBed '{}'", source))?
            .cached();
        let iter = reader.get_interval(chrom, start, end).with_context(|| {
            format!(
                "failed to query BigBed interval {}:{}-{}",
                chrom, start, end
            )
        })?;
        collect_bigbed_features(iter, chrom, limit, &mut features)?;
    } else {
        let path = expand_path(source)?;
        let mut reader = BigBedRead::open_file(&path)
            .with_context(|| format!("failed to open BigBed '{}'", path.display()))?
            .cached();
        let iter = reader.get_interval(chrom, start, end).with_context(|| {
            format!(
                "failed to query BigBed interval {}:{}-{}",
                chrom, start, end
            )
        })?;
        collect_bigbed_features(iter, chrom, limit, &mut features)?;
    }

    Ok(features)
}

fn accumulate_bigwig_records(
    iter: impl Iterator<Item = Result<bigtools::Value, bigtools::BBIReadError>>,
    aggregates: &mut [AggregateBin],
    start: u32,
    end: u32,
    bin_width: u64,
    bins: usize,
) -> Result<()> {
    for record in iter {
        let record = record.context("failed to decode BigWig record")?;
        let overlap_start = u64::from(record.start.max(start));
        let overlap_end = u64::from(record.end.min(end));
        if overlap_end <= overlap_start {
            continue;
        }

        let first_bin = ((overlap_start - u64::from(start)) / bin_width) as usize;
        let last_bin = ((overlap_end - 1 - u64::from(start)) / bin_width) as usize;
        for bin_index in first_bin..=last_bin.min(bins - 1) {
            let bin_start = u64::from(start) + (bin_index as u64 * bin_width);
            let bin_end = (bin_start + bin_width).min(u64::from(end));
            let covered_start = overlap_start.max(bin_start);
            let covered_end = overlap_end.min(bin_end);
            if covered_end <= covered_start {
                continue;
            }
            let overlap = covered_end - covered_start;
            let aggregate = &mut aggregates[bin_index];
            aggregate.weighted_sum += f64::from(record.value) * overlap as f64;
            aggregate.covered_bases += overlap;
            aggregate.max_value = aggregate.max_value.max(f64::from(record.value));
            aggregate.min_value = aggregate.min_value.min(f64::from(record.value));
            aggregate.observed_count += 1;
        }
    }

    Ok(())
}

fn collect_bigbed_features(
    iter: impl Iterator<Item = Result<BedEntry, bigtools::BBIReadError>>,
    chrom: &str,
    limit: usize,
    features: &mut Vec<Feature>,
) -> Result<()> {
    for record in iter.take(limit) {
        let record = record.context("failed to decode BigBed record")?;
        features.push(bigbed_to_feature(chrom, record));
    }
    Ok(())
}

fn bigbed_to_feature(chrom: &str, entry: BedEntry) -> Feature {
    let mut extra_fields = entry.rest.split_whitespace();
    let label = extra_fields.next().map(ToOwned::to_owned);
    let score = extra_fields
        .next()
        .and_then(|value| value.parse::<f32>().ok());
    let strand = extra_fields.next().map(ToOwned::to_owned);
    let attributes = if entry.rest.is_empty() {
        None
    } else {
        Some(entry.rest)
    };

    Feature {
        chrom: chrom.to_string(),
        start: u64::from(entry.start),
        end: u64::from(entry.end),
        label,
        strand,
        score,
        feature_type: Some("bigBed".to_string()),
        attributes,
    }
}

fn parse_bed_line(line: &str) -> Option<Feature> {
    if skip_line(line) {
        return None;
    }

    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 3 {
        return None;
    }

    let start = fields[1].parse::<u64>().ok()?;
    let end = fields[2].parse::<u64>().ok()?;
    let label = fields.get(3).map(|value| (*value).to_string());
    let score = fields.get(4).and_then(|value| value.parse::<f32>().ok());
    let strand = fields.get(5).map(|value| (*value).to_string());
    let attributes = if fields.len() > 6 {
        Some(fields[6..].join("\t"))
    } else {
        None
    };

    Some(Feature {
        chrom: fields[0].to_string(),
        start,
        end,
        label,
        strand,
        score,
        feature_type: Some("bed".to_string()),
        attributes,
    })
}

fn parse_gtf_line(line: &str) -> Option<Feature> {
    if skip_line(line) {
        return None;
    }

    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        return None;
    }

    let start = fields[3].parse::<u64>().ok()?.saturating_sub(1);
    let end = fields[4].parse::<u64>().ok()?;
    let score = if fields[5] == "." {
        None
    } else {
        fields[5].parse::<f32>().ok()
    };
    let strand = if fields[6] == "." {
        None
    } else {
        Some(fields[6].to_string())
    };
    let label = parse_gtf_label(fields[8]);

    Some(Feature {
        chrom: fields[0].to_string(),
        start,
        end,
        label,
        strand,
        score,
        feature_type: Some(fields[2].to_string()),
        attributes: Some(fields[8].trim().to_string()),
    })
}

fn parse_gtf_label(attributes: &str) -> Option<String> {
    for key in ["gene_name", "gene_id", "transcript_name", "transcript_id"] {
        if let Some(value) = parse_gtf_attribute(attributes, key) {
            return Some(value);
        }
    }
    None
}

fn parse_gtf_attribute(attributes: &str, key: &str) -> Option<String> {
    for entry in attributes.split(';') {
        let trimmed = entry.trim();
        if trimmed.is_empty() {
            continue;
        }
        let (name, value) = trimmed.split_once(' ')?;
        if name == key {
            return Some(value.trim_matches('"').to_string());
        }
    }
    None
}

fn skip_line(line: &str) -> bool {
    let trimmed = line.trim();
    trimmed.is_empty()
        || trimmed.starts_with('#')
        || trimmed.starts_with("track ")
        || trimmed.starts_with("browser ")
}

#[cfg(test)]
mod tests {
    use super::{igvjs_format, infer_track_kind, is_genomic_file, parse_bed_line, parse_gtf_line};
    use crate::config::TrackKind;

    #[test]
    fn parse_bed_feature() {
        let feature = parse_bed_line("chr1\t10\t20\tpeak1\t50\t+\textra").unwrap();
        assert_eq!(feature.chrom, "chr1");
        assert_eq!(feature.start, 10);
        assert_eq!(feature.end, 20);
        assert_eq!(feature.label.as_deref(), Some("peak1"));
        assert_eq!(feature.strand.as_deref(), Some("+"));
        assert_eq!(feature.score, Some(50.0));
    }

    #[test]
    fn parse_gtf_feature() {
        let line = "chr1\tsource\texon\t11\t20\t.\t-\t.\tgene_id \"g1\"; gene_name \"GENE1\";";
        let feature = parse_gtf_line(line).unwrap();
        assert_eq!(feature.chrom, "chr1");
        assert_eq!(feature.start, 10);
        assert_eq!(feature.end, 20);
        assert_eq!(feature.label.as_deref(), Some("GENE1"));
        assert_eq!(feature.feature_type.as_deref(), Some("exon"));
        assert_eq!(feature.strand.as_deref(), Some("-"));
    }

    #[test]
    fn is_genomic_file_recognizes_formats() {
        assert!(is_genomic_file("/tmp/sample.bam"));
        assert!(is_genomic_file("/tmp/sample.bam.bai"));
        assert!(is_genomic_file("/tmp/sample.cram"));
        assert!(is_genomic_file("/tmp/sample.vcf.gz"));
        assert!(is_genomic_file("/tmp/sample.bw"));
        assert!(is_genomic_file("/tmp/sample.gff3.gz"));
        assert!(!is_genomic_file("/tmp/readme.txt"));
        assert!(!is_genomic_file("/tmp/data.csv"));
    }

    #[test]
    fn igvjs_format_maps_extensions() {
        assert_eq!(igvjs_format("/tmp/sample.bam"), Some("bam"));
        assert_eq!(igvjs_format("/tmp/sample.bw"), Some("bigwig"));
        assert_eq!(igvjs_format("/tmp/sample.bigWig"), Some("bigwig"));
        assert_eq!(igvjs_format("/tmp/sample.vcf.gz"), Some("vcf"));
        assert_eq!(igvjs_format("/tmp/sample.gff3.gz"), Some("gff3"));
        assert_eq!(igvjs_format("/tmp/sample.txt"), None);
    }

    #[test]
    fn infer_track_kind_by_extension() {
        assert!(matches!(
            infer_track_kind("/tmp/test.bigWig"),
            Some(TrackKind::BigWig)
        ));
        assert!(matches!(
            infer_track_kind("/tmp/test.bb"),
            Some(TrackKind::BigBed)
        ));
        assert!(matches!(
            infer_track_kind("/tmp/test.gtf.gz"),
            Some(TrackKind::Gtf)
        ));
        assert!(matches!(
            infer_track_kind("/tmp/test.bed"),
            Some(TrackKind::Bed)
        ));
    }
}
