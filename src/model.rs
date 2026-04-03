use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::config::{Locus, TrackConfig, TrackKind};

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum WindowFunction {
    None,
    Mean,
    Min,
    Max,
    Count,
    Density,
}

impl Default for WindowFunction {
    fn default() -> Self {
        Self::Mean
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Feature {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub label: Option<String>,
    pub strand: Option<String>,
    pub score: Option<f32>,
    pub feature_type: Option<String>,
    pub attributes: Option<String>,
}

#[derive(Clone, Debug, Serialize)]
pub struct SignalBin {
    pub start: u64,
    pub end: u64,
    pub value: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct TrackResponse {
    pub track_id: String,
    pub kind: TrackKind,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub window_function: Option<WindowFunction>,
    pub signal: Option<Vec<SignalBin>>,
    pub features: Option<Vec<Feature>>,
}

#[derive(Clone, Debug, Serialize)]
pub struct ApiConfigResponse {
    pub title: String,
    pub genome_name: String,
    pub chrom_sizes: BTreeMap<String, u64>,
    pub default_locus: Option<Locus>,
    pub ui: UiConfigResponse,
    pub tracks: Vec<TrackConfig>,
    pub todos: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
pub struct UiConfigResponse {
    pub allowed_roots: Vec<String>,
    pub supports_track_add: bool,
}

#[derive(Clone, Debug, Deserialize)]
pub struct AddTrackRequest {
    pub id: Option<String>,
    pub name: Option<String>,
    pub kind: Option<TrackKind>,
    pub source: String,
    pub style: Option<crate::config::TrackStyle>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct ReorderTracksRequest {
    pub track_ids: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
pub struct FileBrowserEntry {
    pub name: String,
    pub path: String,
    pub is_dir: bool,
    pub size: Option<u64>,
    pub format: Option<String>,
}

#[derive(Clone, Debug, Serialize)]
pub struct FileBrowserResponse {
    pub roots: Vec<String>,
    pub current_path: Option<String>,
    pub parent_path: Option<String>,
    pub entries: Vec<FileBrowserEntry>,
}
