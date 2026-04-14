# genome_viewer

`genome_viewer` is a single-binary Rust web server for browsing genomic tracks in the browser with [igv.js](https://github.com/igvteam/igv.js). It is built for the common bioinformatics case where the data lives on a workstation, lab server, or shared filesystem and you want a browser UI without setting up a separate web stack.

It serves a small embedded frontend, queries BigWig and BigBed server-side, preloads BED and GTF text tracks, and exposes a constrained file browser for loading local tracks at runtime.

## Why this exists

[igv-webapp](https://github.com/igvteam/igv-webapp) is useful for public data and CORS-friendly URLs, but it is still a client-side web app. That means the browser needs direct access to the data source and has no notion of your server filesystem.

`genome_viewer` takes a different approach:

- Single binary, no Node.js, no Java, no separate static asset deployment.
- Zero-config startup for common genomes such as `hg38` and `mm10`.
- Server-side BigWig/BigBed range queries, so the browser receives compact JSON instead of parsing binary formats itself.
- Built-in file browser for loading tracks from allowed local directories.
- Token authentication enabled by default.
- Publication-oriented figure export directly from the browser.

## Highlights

- **Embedded frontend**: `static/index.html` is compiled into the binary, so deployment is just the executable.
- **Layered configuration**: CLI flags, optional JSON config, optional user config, and UCSC chromosome-size fallback.
- **Track management at runtime**: add, remove, and reorder server-managed tracks through the API.
- **Practical browser workflow**: load server files, load remote URLs, save/load sessions, save SVG/PNG, and preserve session state across refresh and re-authentication.
- **Security-minded local file access**: local paths are restricted to canonicalized `allowed_roots`; remote URLs are supported read-only.

## Quick start

### Install

```bash
cargo build --release
# binary: target/release/genome_viewer

# or install into Cargo's bin dir
cargo install --path .
```

### Run

```bash
genome_viewer
genome_viewer --genome mm10
genome_viewer --root ~/data/tracks
genome_viewer --config example-config.json
genome_viewer --bind 0.0.0.0:9000
genome_viewer -p 9000
genome_viewer --port 52000-53000
genome_viewer --token MY_SECRET
genome_viewer --refresh-token
genome_viewer --no-token
```

Default behavior:

- Binds to `0.0.0.0` and picks a random free port in `50000-60000`
- Uses genome `hg38`
- Enables token auth unless `--no-token` is set
- Includes the current working directory as an allowed root unless `--no-cwd` is set

At startup, the server prints a localhost URL plus detected non-loopback network URLs. Open one in your browser and log in with the printed token if auth is enabled.

### Try the bundled example

The repository includes a small demo config and example text tracks:

```bash
cargo run -- --config example-config.json
```

That config demonstrates:

- `genome.default_locus`
- `ui.allowed_roots`
- local BED/GTF tracks
- remote BigWig/BigBed sources

## CLI reference

| Flag | Description |
|------|-------------|
| `--config <path>` | JSON config file. Supports a full config or a tracks-only config. |
| `--genome <name>` | Genome name when not fixed by the JSON config. Default: `hg38`. |
| `--chrom-sizes <path-or-url>` | Chromosome sizes source. Local path or URL. |
| `--bind <addr>` | Bind address. If the bind port is `0` and `--port` is unset, the server chooses a random free port in `50000-60000`. |
| `-p`, `--port <port-or-range>` | Override only the port portion of `--bind`. Accepts either a fixed port like `9000` or an inclusive range like `52000-53000`. |
| `--token [value]` | Enable auth with an explicit token, or auto-generate one if no value is given. |
| `--no-token` | Disable authentication. |
| `--refresh-token` | Generate a new token and save it to `~/.config/genome_viewer/config.yaml`. |
| `--root <path>` | Add an allowed local root. Repeatable. |
| `--no-cwd` | Do not add the current directory to allowed roots. |
| `--title <text>` | Viewer title. |

## Configuration

### User config

Optional defaults live in `~/.config/genome_viewer/config.yaml`:

```yaml
genome: hg38
chrom_sizes: ~/db/gencode/GRCh38/GRCh38.primary_assembly.genome.fa.chromsize
token: abc123def456...
allowed_roots:
  - ~/data/tracks
  - /shared/genomics
```

Use `--refresh-token` to generate and save a new token there. The file is written with `0600` permissions.

### JSON config

Tracks-only config:

```json
{
  "tracks": [
    {
      "id": "signal",
      "name": "My Signal",
      "kind": "bigwig",
      "source": "/path/to/file.bw"
    }
  ]
}
```

Full config:

```json
{
  "title": "My Viewer",
  "genome": {
    "name": "hg38",
    "chrom_sizes": "/path/to/hg38.chrom.sizes",
    "default_locus": {
      "chrom": "chr1",
      "start": 155184000,
      "end": 155194000
    }
  },
  "ui": {
    "allowed_roots": [
      "/data/tracks"
    ]
  },
  "tracks": [
    {
      "id": "peaks",
      "name": "CTCF Peaks",
      "kind": "bed",
      "source": "/path/to/peaks.bed",
      "style": {
        "color": "#2196F3",
        "height": 50
      }
    }
  ]
}
```

Multiple genomes:

```json
{
  "title": "Multi-genome Viewer",
  "genomes": [
    {
      "name": "hg38",
      "label": "Human (hg38)",
      "chrom_sizes": "~/db/gencode/GRCh38/GRCh38.primary_assembly.genome.fa.chromsize"
    },
    {
      "name": "custom_asm",
      "label": "Custom Assembly",
      "chrom_sizes": "/data/custom/custom_asm.chrom.sizes",
      "default_locus": {
        "chrom": "chr1",
        "start": 0,
        "end": 50000
      },
      "reference": {
        "fasta": "/data/custom/custom_asm.fa.gz",
        "fasta_index": "/data/custom/custom_asm.fa.gz.fai",
        "compressed_fasta_index": "/data/custom/custom_asm.fa.gz.gzi",
        "cytoband": "/data/custom/custom_asm.cytoBand.txt",
        "alias": "/data/custom/custom_asm.alias.tsv"
      }
    }
  ],
  "tracks": []
}
```

### Effective precedence

The config is resolved per field, not by one blanket merge order:

- If `--config` contains a `genome` section, that genome block wins over `--genome` and `--chrom-sizes`.
- If `--config` contains `genomes`, those genomes define the genome menu. `--genome` selects the initial genome if it matches one of them.
- Otherwise `--genome` and `--chrom-sizes` override the user config.
- If chromosome sizes are still missing, the server falls back to the UCSC URL for the selected genome.
- If `--config` contains a `title`, it wins over `--title`.
- Allowed roots are merged from:
  - current working directory unless `--no-cwd`
  - repeatable `--root` flags
  - user config `allowed_roots`
  - JSON config `ui.allowed_roots`

## Supported formats

### Server-side queried tracks

These are exposed through `/api/tracks/{id}/query`:

| Format | Extensions | Behavior |
|--------|------------|----------|
| BigWig | `.bw`, `.bigwig` | On-demand query with binning and window functions |
| BigBed | `.bb`, `.bigbed` | On-demand indexed feature query |
| BED | `.bed`, `.bed.gz` | Preloaded into memory at startup or add time |
| GTF | `.gtf`, `.gtf.gz` | Preloaded into memory with 1-based to 0-based conversion |

Supported BigWig window functions:

- `mean`
- `min`
- `max`
- `count`
- `density`
- `none`

### File browser / raw igv.js loading

The built-in file browser also exposes these igv.js-compatible formats as raw files via `/api/data`:

- BAM
- CRAM
- VCF / VCF.GZ
- BED / BED.GZ
- GTF / GTF.GZ
- GFF / GFF3
- WIG
- BedGraph
- SEG
- common index files such as `.bai`, `.crai`, `.tbi`, `.csi`, `.idx`

## Browser workflow

- **Tracks > Server Files...**: browse configured local roots and load genomic files.
- **Tracks > Load URL...**: load a remote igv.js-compatible track by URL.
- **Session > Save Session / Load Session...**: export and restore IGV session JSON.
- **Save Image > Save as SVG / Save as PNG**: use igv.js export for the current browser view.
- **Save Image > Publication Figure...**: generate a cleaner figure-oriented SVG/PNG export from queried track data.

The browser also stores session JSON in `sessionStorage` before refresh or auth redirect, then restores it on the next page load.

## Authentication and local path security

Authentication is enabled by default.

- The server uses a token supplied via `--token`, loaded from user config, or auto-generated at startup.
- Requests are authenticated via the `genome_viewer_token` cookie or `Authorization: Bearer <token>`.
- Unauthenticated requests are redirected to a login page served by the app.
- The auth cookie has `HttpOnly`, `SameSite=Strict`, and `Max-Age=86400`.

Local file access is restricted:

- local paths must resolve inside canonicalized `allowed_roots`
- symlink traversal outside those roots is blocked by canonicalization before validation
- if there are no allowed roots, UI-based local file loading is disabled

Remote `http://` and `https://` track sources are allowed and bypass local path checks, but they are read-only sources.

## API

All endpoints live under `/api/`.

| Method | Path | Description |
|--------|------|-------------|
| `POST` | `/api/auth` | Token login form target |
| `GET` | `/api/config` | Viewer config, chromosome sizes, tracks, UI metadata |
| `GET` | `/api/files?path=` | Browse allowed local directories |
| `GET` | `/api/data?path=` | Serve a local file inside allowed roots |
| `POST` | `/api/tracks` | Add a runtime track |
| `DELETE` | `/api/tracks/{id}` | Remove a runtime track |
| `POST` | `/api/tracks/reorder` | Reorder all server-managed tracks |
| `GET` | `/api/tracks/{id}/query?...` | Query signal or features for a region |

### Add a track

```bash
curl -X POST http://localhost:<port>/api/tracks \
  -H 'Content-Type: application/json' \
  -H 'Authorization: Bearer YOUR_TOKEN' \
  -d '{
    "source": "/data/tracks/sample.bigWig",
    "name": "Sample signal",
    "kind": "bigwig"
  }'
```

`kind` is optional if it can be inferred from the source extension.

### Query a track

```bash
curl "http://localhost:<port>/api/tracks/sample-bigwig/query?chrom=chr1&start=100000&end=120000&bins=800&window_function=mean" \
  -H 'Authorization: Bearer YOUR_TOKEN'
```

Query parameters:

- `chrom` required
- `start` required, 0-based
- `end` required, 0-based exclusive
- `bins` optional, default `800` for BigWig
- `limit` optional, default `2000` for feature tracks
- `window_function` optional, default `mean`

## Publication figure export

The publication figure workflow renders a separate SVG, then optionally rasterizes it to PNG. It is intended for cleaner static figures than the raw igv.js screenshot.

Current figure support includes:

- BigWig signal area plots
- BED / BigBed interval tracks
- GTF gene models with exons, UTR/CDS structure, intron lines, strand arrows, and labels
- coordinate axis
- scale bar
- region label
- ideogram

Configurable figure settings include width in mm, track heights, margins, spacing, font family, font size, label placement, per-track colors, and output DPI.

The figure generator can use:

- server-managed tracks from the backend
- tracks loaded from the server file browser after auto-registering them through `POST /api/tracks`
- some igv.js-managed feature tracks by reading their viewport cache

## Chromosome sizes

Chromosome sizes are resolved as follows:

1. If `--config` contains `genome.chrom_sizes`, use that.
2. Otherwise `--chrom-sizes` overrides the user config value.
3. Otherwise use `chrom_sizes` from `~/.config/genome_viewer/config.yaml`.
4. Otherwise fall back to the UCSC chrom sizes URL for the selected genome.

For common assemblies such as `hg38`, `hg19`, `mm10`, and `mm39`, zero-config startup usually works out of the box.

## Development

Build:

```bash
cargo build
cargo build --release
```

Test:

```bash
cargo test
cargo test slugify
cargo test -- --nocapture
```

Logging:

```bash
RUST_LOG=genome_viewer=debug cargo run
RUST_LOG=genome_viewer=debug,tower_http=debug cargo run
```

## Architecture snapshot

The project is intentionally small:

- `src/main.rs`: CLI, Axum router, handlers, auth middleware, runtime state
- `src/config.rs`: config loading, path normalization, chromosome-size loading, safe source reading
- `src/model.rs`: API request/response types
- `src/tracks.rs`: format inference, text-track parsing, BigWig/BigBed query logic
- `static/index.html`: single-file frontend UI

Key implementation details:

- BigWig and BigBed access is synchronous in `bigtools`, so queries run in `spawn_blocking`.
- BED and GTF tracks are parsed into in-memory per-chromosome vectors and queried by binary-search-like partitioning.
- Internal server errors are logged server-side and returned as a generic `internal server error` message to the client.
- The app prints clickable local and network URLs at startup.
