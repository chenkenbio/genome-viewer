use std::collections::HashMap;

use anyhow::{Context, Result, anyhow, bail};
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

/// Error marker for conditions caused by the *request/data*, not the server —
/// e.g. an HDF5 file in an unsupported layout, a remote HDF5 URL, or a missing
/// chromosome. Carried inside `anyhow::Error` and downcast in `main.rs` so the
/// client gets an HTTP 400 with the explanation instead of a generic 500.
#[derive(Debug)]
pub struct BadRequest(pub String);

impl std::fmt::Display for BadRequest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for BadRequest {}

/// Chromosome-name candidates to probe in an HDF5 file, in priority order:
/// the requested name first, then the name with its `chr` prefix toggled
/// (`chr1` ↔ `1`). seedat `BigWigH5` files name one dataset per chromosome, but
/// the prefix convention varies between sources, so we try both.
fn h5_chrom_candidates(chrom: &str) -> [String; 2] {
    let toggled = match chrom.strip_prefix("chr") {
        Some(rest) => rest.to_string(),
        None => format!("chr{chrom}"),
    };
    [chrom.to_string(), toggled]
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
    } else if lower.ends_with(".h5") || lower.ends_with(".hdf5") {
        Some(TrackKind::Hdf5)
    } else {
        None
    }
}

/// igv.js-recognized track formats, mapping a lowercase filename suffix to the
/// igv.js format name. The format is sent explicitly when a file is loaded,
/// because our data URLs carry the path in a query string (`/api/data?path=…`)
/// that igv.js can't parse for an extension. This is the single source of truth
/// shared by [`igvjs_format`] (the format hint) and [`is_genomic_file`] (browser
/// visibility), so the two never drift apart. Matching is first-wins over this
/// order, but the suffixes are mutually exclusive under `ends_with`, so order is
/// not load-bearing.
const IGV_FORMATS: &[(&str, &str)] = &[
    // Alignments
    (".bam", "bam"),
    (".cram", "cram"),
    // Variants
    (".vcf.gz", "vcf"),
    (".vcf", "vcf"),
    // Signal
    (".bw", "bigwig"),
    (".bigwig", "bigwig"),
    (".bedgraph.gz", "bedgraph"),
    (".bedgraph", "bedgraph"),
    (".bdg", "bedgraph"),
    (".wig", "wig"),
    (".tdf", "tdf"),
    // Indexed binary features
    (".bb", "bigbed"),
    (".bigbed", "bigbed"),
    // Text features
    (".bed.gz", "bed"),
    (".bed", "bed"),
    (".gtf.gz", "gtf"),
    (".gtf", "gtf"),
    (".gff3.gz", "gff3"),
    (".gff3", "gff3"),
    (".gff.gz", "gff"),
    (".gff", "gff"),
    // ENCODE peak formats (BED-derived)
    (".narrowpeak.gz", "narrowpeak"),
    (".narrowpeak", "narrowpeak"),
    (".broadpeak.gz", "broadpeak"),
    (".broadpeak", "broadpeak"),
    (".gappedpeak.gz", "gappedpeak"),
    (".gappedpeak", "gappedpeak"),
    (".regionpeak.gz", "regionpeak"),
    (".regionpeak", "regionpeak"),
    // Interactions
    (".bedpe.gz", "bedpe"),
    (".bedpe", "bedpe"),
    (".interact", "interact"),
    // Gene-prediction tables
    (".refgene", "refgene"),
    (".genepredext", "genepredext"),
    (".genepred", "genepred"),
    (".refflat", "refflat"),
    // Copy-number segments
    (".seg.gz", "seg"),
    (".seg", "seg"),
    // GWAS / mutation
    (".gwas", "gwas"),
    (".maf", "maf"),
    (".mut", "mut"),
    // Methylation
    (".bedmethyl", "bedmethyl"),
];

/// Non-igv-native suffixes the browser should still surface: HDF5 signal
/// (served via a custom server-side track) and sidecar index files.
const EXTRA_VISIBLE_SUFFIXES: &[&str] =
    &[".h5", ".hdf5", ".bai", ".crai", ".tbi", ".csi", ".idx"];

/// Maps a filename to the igv.js format name, or `None` if unrecognized.
pub fn igvjs_format(source: &str) -> Option<&'static str> {
    let lower = source.to_ascii_lowercase();
    IGV_FORMATS
        .iter()
        .find(|(suffix, _)| lower.ends_with(suffix))
        .map(|(_, format)| *format)
}

/// Returns true for any file the browser should show: every igv.js-supported
/// format plus HDF5 and sidecar index files.
pub fn is_genomic_file(source: &str) -> bool {
    let lower = source.to_ascii_lowercase();
    igvjs_format(&lower).is_some() || EXTRA_VISIBLE_SUFFIXES.iter().any(|s| lower.ends_with(s))
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

    let width = u64::from(end - start);
    let bins = effective_signal_bins(bins, width);
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

    Ok(finalize_bins(
        aggregates,
        u64::from(start),
        u64::from(end),
        bin_width,
        window_function,
    ))
}

/// Return the number of signal bins to emit for a query span.
///
/// The requested bin count is capped for payload size, but it also must not
/// exceed the number of bases in the requested span. Otherwise a small zoom
/// window like 100 bp with `bins=4000` produces 100 real 1-bp bins followed by
/// thousands of zero-width/out-of-window bins that clients can mis-render.
fn effective_signal_bins(requested_bins: usize, width: u64) -> usize {
    debug_assert!(width > 0);
    let requested_bins = requested_bins.clamp(1, 4000);
    requested_bins.min(width.min(4000) as usize)
}

/// Turn accumulated per-bin statistics into `SignalBin`s, applying the chosen
/// window function. Shared by `query_bigwig` and `query_hdf5` so both formats
/// produce identical bin semantics (empty bins → 0; mean over covered bases;
/// min/max over observed values; density = covered fraction of the bin span).
fn finalize_bins(
    aggregates: Vec<AggregateBin>,
    start: u64,
    end: u64,
    bin_width: u64,
    window_function: WindowFunction,
) -> Vec<SignalBin> {
    aggregates
        .into_iter()
        .enumerate()
        .filter_map(|(index, aggregate)| {
            let offset = (index as u64).saturating_mul(bin_width);
            let bin_start = start.saturating_add(offset);
            if bin_start >= end {
                return None;
            }
            let bin_end = bin_start.saturating_add(bin_width).min(end);
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
            Some(SignalBin {
                start: bin_start,
                end: bin_end,
                value,
            })
        })
        .collect()
}

/// Query a seedat `BigWigH5` base-resolution signal file.
///
/// Layout (verified on disk): one 1-D HDF5 dataset per chromosome at the file
/// root, dataset name == chromosome name, length == chromosome size in bp,
/// index `i` == 0-based genomic position. Stored dtype is typically `float16`
/// (also handles `float32`/`float64`). Absent positions are stored as NaN.
///
/// The requested window is read in bounded slabs and folded into the same
/// per-bin aggregates as [`query_bigwig`], so both formats share identical
/// binning/window-function semantics. NaN is treated as *no coverage* (skipped),
/// which makes empty/partial bins behave exactly like BigWig absent intervals.
pub fn query_hdf5(
    source: &str,
    chrom: &str,
    start: u64,
    end: u64,
    bins: usize,
    window_function: WindowFunction,
) -> Result<Vec<SignalBin>> {
    if is_remote_path(source) {
        bail!(BadRequest(
            "HDF5 sources must be local files; HTTP-range HDF5 is not supported".to_string()
        ));
    }
    if end <= start {
        bail!("end must be greater than start");
    }

    let width = end - start;
    let bins = effective_signal_bins(bins, width);
    let bin_width = width.div_ceil(bins as u64).max(1);

    let path = expand_path(source)?;
    let file = hdf5::File::open(&path)
        .with_context(|| format!("failed to open HDF5 file '{}'", path.display()))?;

    // Reject the binned LowResBigWigH5 layout (root `__resolution__` marker)
    // before reading anything, with a message the UI can show.
    if let Ok(members) = file.member_names() {
        if members.iter().any(|name| name == "__resolution__") {
            bail!(BadRequest(format!(
                "HDF5 file '{}' uses the binned LowResBigWigH5 layout, which is not supported; \
                 provide a base-resolution BigWigH5 file",
                source
            )));
        }
    }

    // Resolve the per-chromosome dataset: requested name first, then the name
    // with its `chr` prefix toggled.
    let dataset = h5_chrom_candidates(chrom)
        .iter()
        .find_map(|name| file.dataset(name).ok())
        .ok_or_else(|| {
            anyhow!(BadRequest(format!(
                "chromosome '{}' not found in HDF5 file '{}'",
                chrom, source
            )))
        })?;

    if dataset.ndim() != 1 {
        bail!(BadRequest(format!(
            "HDF5 dataset for chromosome '{}' is {}-D; expected a 1-D base-resolution array",
            chrom,
            dataset.ndim()
        )));
    }
    let n = dataset.shape()[0] as u64;
    let elem_size = dataset.dtype()?.size();

    let mut aggregates = vec![
        AggregateBin {
            min_value: f64::INFINITY,
            ..AggregateBin::default()
        };
        bins
    ];

    // Clamp the read window to the dataset bounds; positions past the dataset
    // end carry no coverage (their bins simply stay empty). The returned bin
    // coordinates stay on the *requested* coordinate system regardless.
    let read_start = start.min(n);
    let read_end = end.min(n);

    if read_end > read_start {
        // Bounded slabs: cap peak allocation (a chr1 f16 dataset is ~480 MB raw,
        // so reading a whole-chromosome window at once could exhaust memory).
        const SLAB: u64 = 4_000_000;
        let mut read_ns: u128 = 0;
        let mut bin_ns: u128 = 0;
        let mut slab_lo = read_start;
        while slab_lo < read_end {
            let slab_hi = (slab_lo + SLAB).min(read_end);
            accumulate_hdf5_slab(
                &dataset, elem_size, slab_lo, slab_hi, start, bin_width, bins, &mut aggregates,
                &mut read_ns, &mut bin_ns,
            )?;
            slab_lo = slab_hi;
        }
        tracing::debug!(
            chrom,
            bases = read_end - read_start,
            read_ms = (read_ns / 1_000_000) as u64,
            bin_ms = (bin_ns / 1_000_000) as u64,
            "hdf5 slab read/bin timing"
        );
    }

    Ok(finalize_bins(aggregates, start, end, bin_width, window_function))
}

/// Read one slab `[slab_lo, slab_hi)` of a chromosome dataset and fold its
/// finite values into `aggregates`, dispatching on the stored element byte
/// size: 2 → `f16`, 4 → `f32`, 8 → `f64`. Accumulates the read vs. bin time
/// (nanoseconds) into `read_ns`/`bin_ns` for diagnostics.
fn accumulate_hdf5_slab(
    dataset: &hdf5::Dataset,
    elem_size: usize,
    slab_lo: u64,
    slab_hi: u64,
    start: u64,
    bin_width: u64,
    bins: usize,
    aggregates: &mut [AggregateBin],
    read_ns: &mut u128,
    bin_ns: &mut u128,
) -> Result<()> {
    let lo = slab_lo as usize;
    let hi = slab_hi as usize;
    let read_at = std::time::Instant::now();
    match elem_size {
        2 => {
            let array = dataset
                .read_slice_1d::<half::f16, _>(lo..hi)
                .with_context(|| format!("failed to read HDF5 slab {}..{}", lo, hi))?;
            *read_ns += read_at.elapsed().as_nanos();
            let data = array.as_slice().context("HDF5 slab is not contiguous")?;
            let bin_at = std::time::Instant::now();
            bin_slab(data, |v| v.to_f64(), slab_lo, start, bin_width, bins, aggregates);
            *bin_ns += bin_at.elapsed().as_nanos();
        }
        4 => {
            let array = dataset
                .read_slice_1d::<f32, _>(lo..hi)
                .with_context(|| format!("failed to read HDF5 slab {}..{}", lo, hi))?;
            *read_ns += read_at.elapsed().as_nanos();
            let data = array.as_slice().context("HDF5 slab is not contiguous")?;
            let bin_at = std::time::Instant::now();
            bin_slab(data, |v| f64::from(v), slab_lo, start, bin_width, bins, aggregates);
            *bin_ns += bin_at.elapsed().as_nanos();
        }
        8 => {
            let array = dataset
                .read_slice_1d::<f64, _>(lo..hi)
                .with_context(|| format!("failed to read HDF5 slab {}..{}", lo, hi))?;
            *read_ns += read_at.elapsed().as_nanos();
            let data = array.as_slice().context("HDF5 slab is not contiguous")?;
            let bin_at = std::time::Instant::now();
            bin_slab(data, |v| v, slab_lo, start, bin_width, bins, aggregates);
            *bin_ns += bin_at.elapsed().as_nanos();
        }
        other => bail!(BadRequest(format!(
            "unsupported HDF5 element size {} bytes (expected 2/4/8 for f16/f32/f64)",
            other
        ))),
    }
    Ok(())
}

/// Fold one contiguous slab of base-resolution values into the per-bin
/// aggregates. `data[i]` is the value at genomic position `slab_lo + i`; `to_f64`
/// converts the stored element type (`f16`/`f32`/`f64`) to `f64`.
///
/// Bins are processed as contiguous sub-slices of the slab, so the only integer
/// division happens once per bin (not once per base — the previous per-base
/// `/ bin_width` dominated whole-chromosome queries). The inner accumulation is
/// branchless: NaN is excluded from sum/count via `is_finite()`, while max/min
/// rely on `f64::max`/`f64::min` ignoring NaN, which lets the compiler
/// autovectorize the scan.
fn bin_slab<T: Copy>(
    data: &[T],
    to_f64: impl Fn(T) -> f64,
    slab_lo: u64,
    start: u64,
    bin_width: u64,
    bins: usize,
    aggregates: &mut [AggregateBin],
) {
    if data.is_empty() {
        return;
    }
    let slab_hi = slab_lo + data.len() as u64;
    let first_bin = ((slab_lo - start) / bin_width) as usize;
    let last_bin = (((slab_hi - 1) - start) / bin_width) as usize;

    for bin in first_bin..=last_bin.min(bins - 1) {
        let bin_lo = (start + bin as u64 * bin_width).max(slab_lo);
        let bin_hi = (start + (bin as u64 + 1) * bin_width).min(slab_hi);
        if bin_hi <= bin_lo {
            continue;
        }
        let from = (bin_lo - slab_lo) as usize;
        let to = (bin_hi - slab_lo) as usize;

        let mut sum = 0.0_f64;
        let mut covered = 0_u64;
        let mut max_value = f64::NEG_INFINITY;
        let mut min_value = f64::INFINITY;
        for &raw in &data[from..to] {
            let value = to_f64(raw);
            let finite = value.is_finite();
            sum += if finite { value } else { 0.0 };
            covered += u64::from(finite);
            // `f64::max`/`min` return the non-NaN operand, so NaN is ignored
            // and all-NaN ranges leave the ±INF identities unchanged.
            max_value = max_value.max(value);
            min_value = min_value.min(value);
        }

        let aggregate = &mut aggregates[bin];
        aggregate.weighted_sum += sum;
        aggregate.covered_bases += covered;
        aggregate.observed_count += covered;
        aggregate.max_value = aggregate.max_value.max(max_value);
        aggregate.min_value = aggregate.min_value.min(min_value);
    }
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
    use super::{
        AggregateBin, BadRequest, bin_slab, effective_signal_bins, finalize_bins,
        h5_chrom_candidates, igvjs_format, infer_track_kind, is_genomic_file, parse_bed_line,
        parse_gtf_line, query_hdf5,
    };
    use crate::config::TrackKind;
    use crate::model::WindowFunction;

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
        assert!(is_genomic_file("/tmp/sample.bw.h5"));
        assert!(is_genomic_file("/tmp/sample.hdf5"));
        // Newly recognized igv.js formats
        assert!(is_genomic_file("/tmp/peaks.narrowPeak"));
        assert!(is_genomic_file("/tmp/peaks.broadPeak.gz"));
        assert!(is_genomic_file("/tmp/loops.bedpe"));
        assert!(is_genomic_file("/tmp/study.gwas"));
        assert!(is_genomic_file("/tmp/genes.refGene"));
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
        // gff3 must not be shadowed by the bare `.gff` suffix
        assert_eq!(igvjs_format("/tmp/sample.gff"), Some("gff"));
        assert_eq!(igvjs_format("/tmp/peaks.narrowPeak"), Some("narrowpeak"));
        assert_eq!(igvjs_format("/tmp/peaks.broadPeak.gz"), Some("broadpeak"));
        assert_eq!(igvjs_format("/tmp/loops.bedpe"), Some("bedpe"));
        assert_eq!(igvjs_format("/tmp/genes.genePredExt"), Some("genepredext"));
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
        assert!(matches!(
            infer_track_kind("/tmp/K562.plus.bw.h5"),
            Some(TrackKind::Hdf5)
        ));
        assert!(matches!(
            infer_track_kind("/tmp/test.hdf5"),
            Some(TrackKind::Hdf5)
        ));
    }

    // ── HDF5 binning logic (Tier A: no HDF5 file needed) ────────────

    fn one_empty_bin() -> Vec<AggregateBin> {
        vec![AggregateBin {
            min_value: f64::INFINITY,
            ..AggregateBin::default()
        }]
    }

    #[test]
    fn effective_signal_bins_does_not_exceed_span() {
        assert_eq!(effective_signal_bins(4000, 3), 3);
        assert_eq!(effective_signal_bins(800, 10), 10);
        assert_eq!(effective_signal_bins(0, 10), 1);
        assert_eq!(effective_signal_bins(4000, 10_000), 4000);
    }

    #[test]
    fn finalize_bins_omits_bins_past_window_end() {
        let aggregates = vec![
            AggregateBin {
                weighted_sum: 1.0,
                covered_bases: 1,
                min_value: 1.0,
                max_value: 1.0,
                observed_count: 1,
            };
            4
        ];
        let bins = finalize_bins(aggregates, 10, 12, 1, WindowFunction::Mean);
        assert_eq!(bins.len(), 2);
        assert_eq!((bins[0].start, bins[0].end, bins[0].value), (10, 11, 1.0));
        assert_eq!((bins[1].start, bins[1].end, bins[1].value), (11, 12, 1.0));
    }

    #[test]
    fn finalize_bins_handles_nan_per_window_function() {
        // One bin spanning [0,4); values 1, NaN, 3, NaN at positions 0..4.
        let mut agg = one_empty_bin();
        bin_slab(
            &[1.0_f64, f64::NAN, 3.0, f64::NAN],
            |v| v,
            0, // slab_lo
            0, // start
            4, // bin_width
            1, // bins
            &mut agg,
        );

        let mean = finalize_bins(agg.clone(), 0, 4, 4, WindowFunction::Mean);
        assert_eq!(mean.len(), 1);
        assert_eq!((mean[0].start, mean[0].end), (0, 4));
        assert!((mean[0].value - 2.0).abs() < 1e-9); // (1+3)/2, NaN ignored
        assert_eq!(finalize_bins(agg.clone(), 0, 4, 4, WindowFunction::Max)[0].value, 3.0);
        assert_eq!(finalize_bins(agg.clone(), 0, 4, 4, WindowFunction::Min)[0].value, 1.0);
        // Count = number of finite (covered) bases, NaN excluded.
        assert_eq!(finalize_bins(agg.clone(), 0, 4, 4, WindowFunction::Count)[0].value, 2.0);
        // Density = covered bases / bin span.
        assert!((finalize_bins(agg, 0, 4, 4, WindowFunction::Density)[0].value - 0.5).abs() < 1e-9);
    }

    #[test]
    fn finalize_bins_all_nan_bin_is_zero() {
        let mut agg = one_empty_bin();
        bin_slab(&[f64::NAN, f64::NAN], |v| v, 0, 0, 2, 1, &mut agg);
        for wf in [
            WindowFunction::Mean,
            WindowFunction::Max,
            WindowFunction::Min,
            WindowFunction::Count,
            WindowFunction::Density,
        ] {
            assert_eq!(finalize_bins(agg.clone(), 0, 2, 2, wf)[0].value, 0.0);
        }
    }

    #[test]
    fn h5_chrom_candidates_toggles_chr_prefix() {
        assert_eq!(h5_chrom_candidates("chr1"), ["chr1".to_string(), "1".to_string()]);
        assert_eq!(h5_chrom_candidates("1"), ["1".to_string(), "chr1".to_string()]);
    }

    // ── HDF5 round-trip (Tier B: generated .h5 fixtures) ────────────

    /// Write one chromosome dataset using the same builder pattern as seedat's
    /// BigWigH5 files, optionally with the `lzf` filter (else gzip/deflate).
    fn write_chrom<T: hdf5::H5Type>(file: &hdf5::File, name: &str, values: &[T], lzf: bool) {
        let builder = file.new_dataset_builder();
        let builder = if lzf { builder.lzf() } else { builder.deflate(4) };
        builder
            .with_data(values)
            .create(name)
            .expect("write hdf5 chrom dataset");
    }

    #[test]
    fn query_hdf5_bins_f32_with_lzf_and_nan() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sig_f32.h5");
        {
            let file = hdf5::File::create(&path).unwrap();
            // chr1 length 8; NaN = no coverage.
            let values: [f32; 8] = [1.0, 2.0, f32::NAN, 4.0, f32::NAN, f32::NAN, 7.0, 8.0];
            write_chrom(&file, "chr1", &values, true);
        }
        let src = path.to_str().unwrap();

        let mean = query_hdf5(src, "chr1", 0, 8, 2, WindowFunction::Mean).unwrap();
        assert_eq!(mean.len(), 2);
        assert_eq!((mean[0].start, mean[0].end), (0, 4));
        assert_eq!((mean[1].start, mean[1].end), (4, 8));
        assert!((mean[0].value - 7.0 / 3.0).abs() < 1e-6); // (1+2+4)/3 covered
        assert!((mean[1].value - 7.5).abs() < 1e-6); // (7+8)/2 covered

        let max = query_hdf5(src, "chr1", 0, 8, 2, WindowFunction::Max).unwrap();
        assert_eq!((max[0].value, max[1].value), (4.0, 8.0));
        let count = query_hdf5(src, "chr1", 0, 8, 2, WindowFunction::Count).unwrap();
        assert_eq!((count[0].value, count[1].value), (3.0, 2.0));
    }

    #[test]
    fn query_hdf5_small_window_returns_one_base_bins() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sig_1bp.h5");
        {
            let file = hdf5::File::create(&path).unwrap();
            let values: [f32; 5] = [10.0, 11.0, 12.0, 13.0, 14.0];
            write_chrom(&file, "chr1", &values, true);
        }

        let mean =
            query_hdf5(path.to_str().unwrap(), "chr1", 1, 4, 4000, WindowFunction::Mean).unwrap();
        assert_eq!(mean.len(), 3);
        assert_eq!((mean[0].start, mean[0].end, mean[0].value), (1, 2, 11.0));
        assert_eq!((mean[1].start, mean[1].end, mean[1].value), (2, 3, 12.0));
        assert_eq!((mean[2].start, mean[2].end, mean[2].value), (3, 4, 13.0));
    }

    #[test]
    fn query_hdf5_reads_float16() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sig_f16.h5");
        {
            let file = hdf5::File::create(&path).unwrap();
            let values = [half::f16::from_f32(1.5), half::f16::from_f32(3.5)];
            write_chrom(&file, "chr1", &values, true);
        }
        let mean = query_hdf5(path.to_str().unwrap(), "chr1", 0, 2, 1, WindowFunction::Mean).unwrap();
        assert!((mean[0].value - 2.5).abs() < 1e-3); // f16 precision
    }

    #[test]
    fn query_hdf5_toggles_chr_prefix_and_reads_gzip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sig_nochr.h5");
        {
            let file = hdf5::File::create(&path).unwrap();
            let values: [f32; 4] = [5.0, 5.0, 5.0, 5.0];
            write_chrom(&file, "1", &values, false); // gzip; dataset named "1" (no chr)
        }
        // Query "chr1" resolves to dataset "1" via the chr-prefix toggle.
        let mean = query_hdf5(path.to_str().unwrap(), "chr1", 0, 4, 1, WindowFunction::Mean).unwrap();
        assert_eq!(mean[0].value, 5.0);
    }

    #[test]
    fn query_hdf5_clamps_out_of_range_tail() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sig_clamp.h5");
        {
            let file = hdf5::File::create(&path).unwrap();
            let values: [f32; 8] = [2.0; 8];
            write_chrom(&file, "chr1", &values, true);
        }
        // Window [0,12) over a length-8 chrom, 3 bins of width 4: the last bin
        // (8..12) lies entirely past the dataset → no coverage → 0.
        let mean = query_hdf5(path.to_str().unwrap(), "chr1", 0, 12, 3, WindowFunction::Mean).unwrap();
        assert_eq!(mean.len(), 3);
        assert_eq!(mean[0].value, 2.0);
        assert_eq!(mean[1].value, 2.0);
        assert_eq!(mean[2].value, 0.0);
        assert_eq!((mean[2].start, mean[2].end), (8, 12));
    }

    #[test]
    fn query_hdf5_rejects_binned_layout() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("lowres.h5");
        {
            let file = hdf5::File::create(&path).unwrap();
            let marker: [i32; 1] = [10];
            file.new_dataset_builder()
                .with_data(&marker)
                .create("__resolution__")
                .unwrap();
            let values: [f32; 4] = [1.0; 4];
            write_chrom(&file, "chr1", &values, false);
        }
        let err = query_hdf5(path.to_str().unwrap(), "chr1", 0, 4, 1, WindowFunction::Mean)
            .unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("LowResBigWigH5") || msg.contains("binned"), "got: {msg}");
        // Must be a BadRequest (→ HTTP 400), not a generic internal error.
        assert!(err.downcast_ref::<BadRequest>().is_some());
    }

    #[test]
    fn query_hdf5_missing_chrom_is_bad_request() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sig_missing.h5");
        {
            let file = hdf5::File::create(&path).unwrap();
            let values: [f32; 4] = [1.0; 4];
            write_chrom(&file, "chr1", &values, true);
        }
        let err = query_hdf5(path.to_str().unwrap(), "chrZ", 0, 4, 1, WindowFunction::Mean)
            .unwrap_err();
        assert!(err.downcast_ref::<BadRequest>().is_some());
        assert!(err.to_string().contains("chrZ"));
    }

    #[test]
    fn query_hdf5_rejects_remote_source() {
        let err = query_hdf5(
            "https://example.com/sig.h5",
            "chr1",
            0,
            100,
            10,
            WindowFunction::Mean,
        )
        .unwrap_err();
        assert!(err.downcast_ref::<BadRequest>().is_some());
        assert!(err.to_string().contains("local"));
    }
}
