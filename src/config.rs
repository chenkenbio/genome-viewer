use std::collections::{BTreeMap, HashMap};
use std::io::Read;
use std::path::{Component, Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct AppConfig {
    pub title: String,
    pub genome: GenomeConfig,
    #[serde(default)]
    pub ui: UiConfig,
    pub tracks: Vec<TrackConfig>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct GenomeConfig {
    pub name: String,
    pub chrom_sizes: String,
    pub default_locus: Option<Locus>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Locus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct UiConfig {
    #[serde(default)]
    pub allowed_roots: Vec<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum TrackKind {
    Bed,
    BigBed,
    BigWig,
    Gtf,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct TrackStyle {
    pub color: Option<String>,
    pub height: Option<u32>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct TrackConfig {
    pub id: String,
    pub name: String,
    pub kind: TrackKind,
    pub source: String,
    pub style: Option<TrackStyle>,
}

#[derive(Clone, Debug)]
pub struct ResolvedConfig {
    pub raw: AppConfig,
    pub chrom_sizes: BTreeMap<String, u64>,
    pub track_map: HashMap<String, TrackConfig>,
}

/// Input config format for JSON file. All fields except tracks are optional,
/// allowing a tracks-only config: `{ "tracks": [...] }`
#[derive(Clone, Debug, Deserialize)]
pub struct InputConfig {
    pub title: Option<String>,
    pub genome: Option<GenomeConfig>,
    #[serde(default)]
    pub ui: UiConfig,
    #[serde(default)]
    pub tracks: Vec<TrackConfig>,
}

/// User-level defaults from `~/.config/genome_viewer/config.yaml`
#[derive(Clone, Debug, Default, Deserialize)]
pub struct UserConfig {
    pub genome: Option<String>,
    pub chrom_sizes: Option<String>,
    pub token: Option<String>,
    #[serde(default)]
    pub allowed_roots: Vec<String>,
}

/// Load user config from `~/.config/genome_viewer/config.yaml`.
/// Never fails — returns defaults if the file is missing or invalid.
pub fn load_user_config() -> UserConfig {
    let config_path = match expand_path("~/.config/genome_viewer/config.yaml") {
        Ok(p) => p,
        Err(_) => return UserConfig::default(),
    };
    if !config_path.exists() {
        return UserConfig::default();
    }
    match std::fs::read_to_string(&config_path) {
        Ok(text) => match serde_yaml::from_str(&text) {
            Ok(c) => {
                tracing::info!(path = %config_path.display(), "loaded user config");
                c
            }
            Err(e) => {
                tracing::warn!(path = %config_path.display(), error = %e, "failed to parse user config");
                UserConfig::default()
            }
        },
        Err(e) => {
            tracing::warn!(path = %config_path.display(), error = %e, "failed to read user config");
            UserConfig::default()
        }
    }
}

/// Return the UCSC-hosted chrom sizes URL for a genome assembly.
pub fn default_chrom_sizes_url(genome: &str) -> String {
    format!("https://hgdownload.cse.ucsc.edu/goldenpath/{genome}/bigZips/{genome}.chrom.sizes")
}

/// Load a JSON config file (supports both full and tracks-only format).
pub fn load_input_config(path: &Path) -> Result<InputConfig> {
    let expanded_path = expand_path(path.to_string_lossy().as_ref())?;
    let config_text = std::fs::read_to_string(&expanded_path)
        .with_context(|| format!("failed to read config {}", expanded_path.display()))?;
    serde_json::from_str(&config_text)
        .with_context(|| format!("failed to parse config {}", expanded_path.display()))
}

/// Build a `ResolvedConfig` from individual components.
pub fn build_config(
    title: String,
    genome_name: String,
    chrom_sizes_source: String,
    default_locus: Option<Locus>,
    tracks: Vec<TrackConfig>,
) -> Result<ResolvedConfig> {
    let chrom_sizes = load_chrom_sizes(&chrom_sizes_source)?;
    if chrom_sizes.is_empty() {
        bail!(
            "chromosome sizes file '{}' is empty",
            chrom_sizes_source
        );
    }
    let mut track_map = HashMap::new();
    for track in &tracks {
        if track_map.insert(track.id.clone(), track.clone()).is_some() {
            bail!("duplicate track id '{}'", track.id);
        }
    }
    let raw = AppConfig {
        title,
        genome: GenomeConfig {
            name: genome_name,
            chrom_sizes: chrom_sizes_source,
            default_locus,
        },
        ui: UiConfig::default(),
        tracks,
    };
    Ok(ResolvedConfig {
        raw,
        chrom_sizes,
        track_map,
    })
}

pub fn load_chrom_sizes(source: &str) -> Result<BTreeMap<String, u64>> {
    let text = read_source_to_string(source)?;
    let mut chrom_sizes = BTreeMap::new();

    for (line_number, line) in text.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let mut fields = trimmed.split_whitespace();
        let chrom = fields
            .next()
            .with_context(|| format!("missing chromosome at line {}", line_number + 1))?;
        let size = fields
            .next()
            .with_context(|| format!("missing size at line {}", line_number + 1))?
            .parse::<u64>()
            .with_context(|| format!("invalid size at line {}", line_number + 1))?;
        chrom_sizes.insert(chrom.to_string(), size);
    }

    Ok(chrom_sizes)
}

pub fn expand_path(path: &str) -> Result<PathBuf> {
    let expanded =
        shellexpand::full(path).with_context(|| format!("failed to expand path '{}'", path))?;
    Ok(PathBuf::from(expanded.as_ref()))
}

pub fn normalize_local_path(path: &str) -> Result<PathBuf> {
    let expanded = expand_path(path)?;
    let base = if expanded.is_absolute() {
        expanded
    } else {
        std::env::current_dir()
            .context("failed to resolve current working directory")?
            .join(expanded)
    };
    Ok(normalize_path_components(&base))
}

pub fn normalize_path_components(path: &Path) -> PathBuf {
    let mut normalized = PathBuf::new();
    for component in path.components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                normalized.pop();
            }
            other => normalized.push(other.as_os_str()),
        }
    }
    normalized
}

pub fn is_remote_path(path: &str) -> bool {
    path.starts_with("http://") || path.starts_with("https://")
}

/// Maximum file size for text track loading (500 MB).
/// Prevents memory exhaustion from very large BED/GTF files or gzip bombs.
const MAX_SOURCE_BYTES: u64 = 500 * 1024 * 1024;

pub fn read_source_bytes(source: &str) -> Result<Vec<u8>> {
    if is_remote_path(source) {
        let response = reqwest::blocking::get(source)
            .with_context(|| format!("failed to fetch remote source '{}'", source))?
            .error_for_status()
            .with_context(|| format!("remote source returned an error '{}'", source))?;
        if let Some(len) = response.content_length() {
            if len > MAX_SOURCE_BYTES {
                bail!(
                    "remote source '{}' is too large ({} bytes, limit {})",
                    source,
                    len,
                    MAX_SOURCE_BYTES
                );
            }
        }
        let bytes = response
            .bytes()
            .with_context(|| format!("failed to read remote source '{}'", source))?;
        if bytes.len() as u64 > MAX_SOURCE_BYTES {
            bail!(
                "remote source '{}' exceeds size limit ({} bytes, limit {})",
                source,
                bytes.len(),
                MAX_SOURCE_BYTES
            );
        }
        Ok(bytes.to_vec())
    } else {
        let path = expand_path(source)?;
        let meta = std::fs::metadata(&path)
            .with_context(|| format!("failed to stat '{}'", path.display()))?;
        if meta.len() > MAX_SOURCE_BYTES {
            bail!(
                "source '{}' is too large ({} bytes, limit {})",
                source,
                meta.len(),
                MAX_SOURCE_BYTES
            );
        }
        let bytes = std::fs::read(&path)
            .with_context(|| format!("failed to read local source '{}'", path.display()))?;
        Ok(bytes)
    }
}

pub fn read_source_to_string(source: &str) -> Result<String> {
    let bytes = read_source_bytes(source)?;
    let mut output = String::new();

    if source.ends_with(".gz") {
        let mut decoder = GzDecoder::new(bytes.as_slice());
        decoder
            .read_to_string(&mut output)
            .with_context(|| format!("failed to decompress '{}'", source))?;
        if output.len() as u64 > MAX_SOURCE_BYTES {
            bail!(
                "decompressed source '{}' exceeds size limit ({} bytes, limit {})",
                source,
                output.len(),
                MAX_SOURCE_BYTES
            );
        }
    } else {
        output = String::from_utf8(bytes)
            .with_context(|| format!("source '{}' is not valid UTF-8", source))?;
    }

    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::{load_chrom_sizes, normalize_local_path};

    #[test]
    fn parse_chrom_sizes_text() {
        let temp = tempfile_content("chr1\t100\nchr2\t200\n");
        let chrom_sizes = load_chrom_sizes(temp.to_str().unwrap()).unwrap();
        assert_eq!(chrom_sizes.get("chr1"), Some(&100));
        assert_eq!(chrom_sizes.get("chr2"), Some(&200));
    }

    #[test]
    fn normalize_relative_local_path() {
        let normalized = normalize_local_path("./foo/../bar/test.bed").unwrap();
        assert!(normalized.ends_with("bar/test.bed"));
    }

    fn tempfile_content(contents: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir();
        let file = dir.join(format!("genome_viewer_test_{}.sizes", std::process::id()));
        std::fs::write(&file, contents).unwrap();
        file
    }
}
