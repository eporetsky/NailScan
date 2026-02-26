# NailScan

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

NailScan uses [nail](https://github.com/TravisWheelerLab/nail) (Nail is an Alignment
Inference tooL) to run fast profile HMM searches against InterPro member
databases (Pfam, NCBIfam, SUPERFAMILY, PANTHER, Gene3D).  Nail uses MMseqs2
for seeding and approximates the HMMER3 Forward/Backward algorithm, giving
large speed-ups over a standard HMMER search.

NailScan is **not** meant to reproduce identical results to InterProScan or to
replace it.  It is intended as a faster alternative for HMM-based protein
annotation where throughput and scalability matter more than exact parity.
NailScan is still under development. Benchmarking Arabidopsis InterProScan with
NailScan resulted in an average 0.95 F1 score for the five member databases
in just over 20 minutes.

> **Disclaimer:** NailScan is not affiliated with nail, InterPro, or any member
> database.  Results are post-processed similarly to InterProScan's filtering
> conventions as closely as possible, but many differences exist.

**References**

- nail: [doi:10.1101/2024.01.27.577580](https://doi.org/10.1101/2024.01.27.577580)
- InterPro: [doi:10.1093/nar/gkac993](https://doi.org/10.1093/nar/gkac993)
- Pfam: [doi:10.1093/nar/gkaa913](https://doi.org/10.1093/nar/gkaa913)

---

## Quick start

### 1. Create the conda environment

```bash
conda env create -f environment.yml
conda activate nailscan
```

### 2. Install nail and HMM databases

```bash
./install.sh
```

This will:
- Clone and build **nail** (requires Rust; only needed at install time).
- Download the five enabled HMM databases from the EBI FTP into `data/`.
- Pre-process the **SUPERFAMILY** HMM file (add COMPO lines required by nail,
  remove single-node models).
- Download **InterPro / GO mapping files** for `--iprlookup` / `--goterms`.

To install only specific databases:

```bash
NAILSCAN_DBS="pfam ncbifam" ./install.sh
```

### 3. Run annotation

**Single FASTA → one merged TSV:**

```bash
./nailscan.single.sh -f sequences.fasta -b myrun
# → results/myrun.tsv
```

**Directory of genome FASTAs → one gzipped TSV per genome:**

```bash
./nailscan.batch.sh -b mybatch -f path/to/fasta_dir/
# → results/mybatch/<genome_id>.tsv.gz
```

Both scripts require **`-b`** (output base name). Run with **`--help`** for
all options.

---

## Supported databases

Five InterPro member databases are enabled by default.  Each uses
database-specific post-processing (`scripts/add_ipr_go.py`) to match
InterProScan's filtering approach as closely as possible.

| Database | InterPro label | Filtering method |
|----------|---------------|-----------------|
| **Pfam** | Pfam | GA domain + sequence thresholds; clan-based overlap deduplication |
| **NCBIfam** | NCBIfam | TC sequence threshold |
| **SUPERFAMILY** | SUPERFAMILY | E-value ≤ 1e-3; global per-protein non-overlapping (35% overlap threshold) |
| **PANTHER** | PANTHER | Best hit per protein; family descriptions from PANTHER classifications |
| **Gene3D** | Gene3D | Bitscore ≥ 10, E-value ≤ 1e-4; global per-protein non-overlapping (20% overlap threshold) |

### Databases not currently enabled

Several databases present in InterProScan are not enabled by default, either
because they require algorithms that cannot be replicated without HMMER's full
output, or because the raw HMM search precision appeared too low to be included.
They can still be run with `--appl <db>`, but results will be less reliable.

| Database | Reason disabled |
|----------|----------------|
| **SFLD** | Requires per-residue active-site matching (`sfld_postprocess`) using HMMER alignment files that nail does not produce. Without this filter, ~20–25% over-prediction. |
| **HAMAP** | Uses Prosite-matrix normalised scores (N_SCORE) for thresholding, not HMM bit scores; cutoffs are not transferable, causing ~10× over-prediction. |
| **PIRSF** | Hierarchy-aware scoring (length, coverage, parent-suppression) is not fully replicable without InterProScan internals; over-predicts broad parent families ~4×. |
| **PIRSR** | PIR Site Rules require residue-level functional site detection, not available from a plain HMM search. |
| **Prosite** | Profiles are in `.prf` format; no standard conversion to HMMER3 profile HMMs. |
| **SMART** | Distributed as HMMER2 format; could be converted but not currently tested. |
| **CDD** | Uses RPS-BLAST, not profile HMMs. |
| **PRINTS** | Fingerprint-based, not HMM-based. |
| **AntiFam** | Included in the download but rarely annotates proteins; enabled with `--appl antifam`. |

See [`docs/FILTERING.md`](docs/FILTERING.md) for detailed filtering logic and
[`docs/INTERPRO_DATABASES.md`](docs/INTERPRO_DATABASES.md) for HMM format
compatibility notes.

---

## Performance

### Accuracy vs InterProScan

Benchmarked against InterProScan 5.77-108.0 on the *Arabidopsis thaliana*
Col-0 proteome (35,386 proteins, 72 threads).

| Database | NailScan predictions | IPR reference | TP | FP | FN | Precision | Recall | F1 |
|----------|---------------------:|-------------:|---:|---:|---:|----------:|-------:|---:|
| **Pfam** | 33,753 | 34,835 | 32,654 | 1,099 | 2,181 | 0.967 | 0.937 | 0.952 |
| **NCBIfam** | 4,075 | 4,080 | 4,004 | 71 | 76 | 0.983 | 0.981 | 0.982 |
| **SUPERFAMILY** | 20,170 | 21,494 | 19,735 | 435 | 1,759 | 0.978 | 0.918 | 0.947 |
| **PANTHER** | 26,170 | 24,627 | 24,178 | 1,992 | 449 | 0.924 | 0.982 | 0.952 |
| **Gene3D** | 23,850 | 25,330 | 23,206 | 644 | 2,124 | 0.973 | 0.916 | 0.944 |
| **Overall** | **108,018** | **110,366** | **103,777** | **4,241** | **6,589** | **0.961** | **0.940** | **0.950** |

### HMM database load time: text vs. pre-serialized binary

Nail can pre-serialize HMM databases to a binary format (`.nailbin`) on first
use.  Subsequent runs load from the binary, which is significantly faster for
large databases.

| Database | Text HMM (s) | Binary cache (s) | Speedup |
|----------|-------------:|-----------------:|--------:|
| Gene3D | 136.1 | 6.7 | 20× |
| PANTHER | 89.4 | 4.4 | 20× |
| NCBIfam | 70.3 | 3.5 | 20× |
| Pfam | 55.9 | 2.7 | 21× |
| SUPERFAMILY | 33.0 | 3.5 | 9× |
| AntiFam | 0.5 | < 0.1 | — |

The binary cache is written automatically on first use; no manual step is
required.

---

## Output format

Both scripts produce a tab-separated TSV with one row per domain match:

| Column | Description |
|--------|-------------|
| `target` | Protein / gene identifier |
| `Analysis` | Database name (e.g. `Pfam`, `Gene3D`, `SUPERFAMILY`) |
| `ACC` | Canonical accession for the matched model |
| `DESC` | Human-readable description |
| `target_start` | Start position on the protein (1-based) |
| `target_end` | End position on the protein (inclusive) |
| `query_start` | Start of the alignment on the HMM |
| `query_end` | End of the alignment on the HMM |
| `score` | Bit score |
| `bias` | Composition bias |
| `evalue` | E-value |
| `cell_frac` | Fraction of DP matrix cells computed by nail |
| `InterPro` | InterPro entry ID *(requires `--iprlookup`)* |
| `IPR_desc` | Short description of the InterPro entry *(requires `--iprlookup`)* |
| `GO` | GO term IDs, semicolon-separated *(requires `--goterms`)* |

Rows are sorted by `target` (gene ID). Single mode writes
`results/<base>.tsv`; batch mode writes
`results/<base>/<genome_id>.tsv.gz`.

---

## Options

```bash
./nailscan.single.sh -f <input.fasta> -b <base> [-o output_dir] [-t threads]
                     [--appl pfam,ncbifam,...] [--iprlookup] [--goterms]

./nailscan.batch.sh  -b <base> [-f fasta_dir] [-o results_dir] [-t threads]
                     [--batch-size N] [--appl pfam,ncbifam,...] [--iprlookup] [--goterms]
```

| Flag | Description | Default |
|------|-------------|---------|
| `-f` | Input FASTA (single) or directory of FASTAs (batch) | `fasta/` |
| `-b` | Output base name **(required)** | — |
| `-o` | Output root directory | `results/` |
| `-t` | Threads | `nproc` |
| `-d` | Data directory | `data/` |
| `--appl` | Comma-separated list of databases to run | all enabled |
| `--batch-size` | Genomes per nail call in batch mode (0 = all) | 0 |
| `--iprlookup` | Add `InterPro` and `IPR_desc` columns | off |
| `--goterms` | Add `GO` column | off |

---

## Configuration

**`config.json`** lists all databases with their EBI FTP URLs, MD5 checksums,
and expected HMM filenames.  Databases with `"enabled": false` are skipped by
default but can still be downloaded and run with `--appl`.

**`install.sh`** reads `config.json` to download and prepare each database.
Run `./install.sh` from the repo root with the `nailscan` conda environment
active.

---

## Repository layout

| Path | Purpose |
|------|---------|
| `config.json` | Database download URLs and enabled/disabled status |
| `install.sh` | Downloads databases; builds nail; sets up mapping files |
| `nailscan.single.sh` | Annotates a single FASTA; writes one merged TSV |
| `nailscan.batch.sh` | Wrapper that calls `nailscan.single.sh` per batch of genomes |
| `environment.yml` | Conda environment (Python, Rust, MMseqs2, HMMER) |
| `scripts/add_ipr_go.py` | Database-specific filtering and InterPro/GO annotation |
| `scripts/benchmark.py` | Compare NailScan output against InterProScan reference |
| `scripts/download_interpro_mappings.sh` | Downloads InterPro/GO mapping files |
| `scripts/cath-resolve-hits` | Gene3D DP domain architecture solver (from InterProScan) |
| `scripts/assign_cath_superfamilies.py` | Gene3D G3DSA mapping script (from InterProScan) |
| `scripts/sfld_postprocess` | SFLD active-site filter binary (from InterProScan; requires HMMER output) |
| `docs/FILTERING.md` | Detailed filtering logic for each database |
| `docs/INTERPRO_DATABASES.md` | HMM format compatibility notes |

---

## Contributors

- **Elly Poretsky ([@eporetsky](https://github.com/eporetsky))**
- Contributions and collaborations welcome - open an issue or pull request.
