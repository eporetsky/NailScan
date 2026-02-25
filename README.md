# NailScan

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

NailScan uses [nail](https://github.com/TravisWheelerLab/nail) (an Alignment Inference tooL) to run fast profile HMM search against HMM-based databases such as [Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/), NCBIFam, Superfamily, and others from InterPro. Nail uses MMseqs2 for seeding and approximates the HMMER3 Forward/Backward algorithm for speed. NailScan is **not** meant to reproduce identical results to InterProScan or to replace it; it is intended as an alternative for HMM-based analyses where speed and scalability are priorities.

[![InterPro](https://img.shields.io/badge/InterPro–Pfam-https%3A%2F%2Fwww.ebi.ac.uk%2Finterpro%2F-blue?labelColor=0a5f9e)](https://www.ebi.ac.uk/interpro/)
[![nail](https://img.shields.io/badge/nail-https%3A%2F%2Fgithub.com%2FTravisWheelerLab%2Fnail-green?labelColor=2ea043)](https://github.com/TravisWheelerLab/nail)

> **Disclaimer:** NailScan is not affiliated with nail, InterPro, or Pfam. This project is independent of the nail developers, the InterPro/Pfam resource maintainers, and the EBI. NailScan is experimental and should be used with caution; results may need to be validated or filtered (e.g. by E‑value) depending on your use case.

**References**

- nail (search engine): [doi:10.1101/2024.01.27.577580](https://doi.org/10.1101/2024.01.27.577580)
- Pfam: [doi:10.1093/nar/gkaa913](https://doi.org/10.1093/nar/gkaa913)
- InterPro: [doi:10.1093/nar/gkac993](https://doi.org/10.1093/nar/gkac993)

---

## Quick start

1. **Create and activate the conda environment**

   ```bash
   conda env create -f environment.yml
   conda activate nailscan
   ```

2. **Install nail and HMM databases**

   ```bash
   ./install.sh
   ```

   With the nailscan env active, this will:
   - Clone and build **nail** (requires **Rust** in the env; only needed at install time).
   - Download HMM databases from InterPro/EBI per `config.json` (Pfam, NCBIFam, Superfamily, AntiFam, HAMAP, SFLD, PIRSF, PIRSR, PANTHER, CATH) into `data/<dbname>/`. Existing DBs are skipped.
   - Convert **Superfamily** from HMMER2 to HMMER3 format if needed (requires **hmmer** in the env).
   - Download **InterPro/GO** mapping files into `data/` if missing (for `--iprlookup` / `--goterms`).

3. **Run search**

   - **Single FASTA** (one output file combining all selected DBs):
     ```bash
     ./nailscan.single.sh -f path/to/sequences.fasta -b myrun
     ```
     Writes **`results/myrun.tsv`** (one merged TSV with an `Analysis` column for the database name).

   - **Batch** (directory of FASTA files):
     ```bash
     ./nailscan.batch.sh -b mybatch -f path/to/fasta_dir
     ```
     Writes **`results/mybatch/<genome_id>.tsv.gz`** (one merged, gzipped TSV per genome).

   Both scripts require **`-b` / `--output-file-base`** to name the output. Use **`--help`** for full options.

### Optional: test run

The `test/` directory contains small FASTA files. After installation:

**Single FASTA:**

```bash
./nailscan.single.sh -f test/test1.fa -b test1 -o results
# → results/test1.tsv (all DBs merged, with Analysis column)
```

**Batch:**

```bash
./nailscan.batch.sh -b testbatch -f test -o results
# → results/testbatch/test1.tsv.gz, results/testbatch/test2.tsv.gz, ...
```

### Python CLI (optional)

Same workflows are available via a Click CLI. Run from the repo root with the nailscan env active and nail on PATH:

```bash
conda activate nailscan
PYTHONPATH=. python -m nailscan single -f test/test1.fa -b test1
PYTHONPATH=. python -m nailscan batch -b mybatch -f test/
```

Use `python -m nailscan single --help` and `python -m nailscan batch --help` for options. This project is not pip-installable; use the nailscan conda env.

## Threading and resources

- Use **`-t N`** to set the number of threads (default: `nproc` or 16).

  ```bash
  ./nailscan.single.sh -f my.faa -b myrun -t 72
  ./nailscan.batch.sh -b mybatch -f /path/to/fastas -t 72
  ```

- **Memory:** Batch mode can use a lot of RAM. For very large runs (e.g. millions of sequences), plan resources accordingly.

## Default behaviour and options

- **By default both scripts run all databases** in `config.json`: pfam, ncbifam, superfamily, antifam, hamap, sfld, pirsf, pirsr, panther, cath. Use **`--applications`** to run only selected DBs, e.g. `--applications pfam,ncbifam`.
- **Single mode** can run one HMM file with **`-h <path>`** (e.g. `-h data/pfam/pfam_a.hmm`); output is still one merged TSV named by `-b`.
- **`--iprlookup`** adds an InterPro column; **`--goterms`** adds a GO terms column. Both need mapping files in `data/` (created by `install.sh` or **`./scripts/download_interpro_mappings.sh`**). See **`docs/INTERPRO_DATABASES.md`** for details.

| Script               | Required options      | Output |
|----------------------|------------------------|--------|
| `nailscan.single.sh` | **`-f`** FASTA, **`-b`** output base | **`<output_dir>/<output_file_base>.tsv`** (one merged TSV, all DBs; `Analysis` column = DB name) |
| `nailscan.batch.sh`  | **`-b`** output base   | **`<output_dir>/<output_file_base>/<genome_id>.tsv.gz`** (one merged TSV per genome) |

Options: **`-f`** input FASTA (single) or directory (batch), **`-b`** / **`--output-file-base`** output name (required), **`-o`** output directory, **`-h`** single HMM path (single only), **`-t`** threads, **`--applications`** comma-separated DB list, **`--iprlookup`**, **`--goterms`**. Run with **`--help`** for full usage.

```bash
./nailscan.single.sh -f <input.fasta> -b <output_base> [-o output_dir] [-t N] [--applications pfam,ncbifam] [--iprlookup] [--goterms]
./nailscan.batch.sh -b <output_base> [-f fasta_dir] [-o results_dir] [-t N] [--applications pfam,ncbifam] [--iprlookup] [--goterms]
```

## Output format

- **nailscan.single.sh** writes one merged TSV: **`<output_dir>/<output_file_base>.tsv`**. Columns include **`Analysis`** (database name), `target`, `NAME`, `ACC`, `DESC`, `target_start`, `target_end`, `query_start`, `query_end`, `score`, `bias`, `evalue`, `cell_frac`, and optionally `InterPro` and `GO` when **`--iprlookup`** / **`--goterms`** are used.

- **nailscan.batch.sh** writes **`<output_dir>/<output_file_base>/<genome_id>.tsv.gz`**: same column layout, one gzipped merged TSV per genome.

**E‑value filtering:** Consider filtering by E‑value (and/or score) in downstream analysis.

### PANTHER output

PANTHER results are automatically post-processed by `scripts/add_ipr_go.py` (regardless of `--iprlookup`/`--goterms`):

1. **ID cleanup** — model names are stripped of the `.orig.30.pir` suffix present in PANTHER's HMM files; `PTHR22842.orig.30.pir` → `PTHR22842`.
2. **Best-hit filter** — InterProScan reports at most one PANTHER family per protein (the highest-scoring family match). NailScan applies the same filter, collapsing typically ~20× more raw hits to one per protein.
3. **Family description** — the `DESC` column is populated from `data/panther/PANTHER*_HMM_classifications`, downloaded by `install.sh`.

Raw PANTHER results are therefore not written to the final TSV; only the best-hit, annotated rows are included.

## Configuration

- **`config.json`** — Used by `install.sh` to download HMM databases from the [EBI iprscan FTP](https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/6/6.0/). Each entry has `url`, `md5_url`, and `hmm_name` (the main HMM file after unpack). Supported DBs: **pfam**, **ncbifam**, **superfamily**, **antifam**, **hamap**, **sfld**, **pirsf**, **pirsr**, **panther**, **cath**. Run `./install.sh` to install all, or e.g. **`NAILSCAN_DBS="pfam antifam" ./install.sh`** to install only some. Superfamily: `install.sh` converts to HMMER3 if needed, removes model_length=1 models, then **replaces** `hmmlib_1.75` with this processed file (original saved as `hmmlib_1.75.raw`). So the config path always points at the nail-safe HMM. See **`docs/INTERPRO_DATABASES.md`** for compatibility notes.

### Model length 1 and `data/install.log`

Some HMM databases (e.g. Superfamily) include **single-node** profile HMMs (**model_length = 1**): the profile has only one match state, so it describes a single alignment column rather than a multi-position motif. They are valid in HMMER3 and **InterProScan does not skip them**; HMMER3 can run them. **Nail** instead assumes at least two nodes and would crash on them. To avoid patching nail, `install.sh` runs **`scripts/filter_length1_hmms.py`** on the Superfamily HMM after conversion: it removes any model with `LENG 1` and **appends a line to `data/install.log`** listing the removed model names. That way the installed HMM file is safe for nail and you have a record of what was dropped. If no models are removed, nothing is appended to the log.

## What each file does

| File | Purpose |
|------|--------|
| `environment.yml` | Conda env: **python**, **rust** (to build nail at install time), **mmseqs2** (runtime), **hmmer** (Superfamily conversion), **click** (optional CLI). |
| `config.json` | DB download URLs and HMM filenames for `install.sh`. |
| `install.sh` | Builds nail, installs it into the env; downloads HMM DBs to `data/<dbname>/` (with MD5 check); for Superfamily, builds a nail-ready HMM (HMMER3 + no model_length=1) and **replaces** `hmmlib_1.75` with it (original → `hmmlib_1.75.raw`), logs removed models to **`data/install.log`**; fetches InterPro/GO mappings if missing. Defines **`make_pfam_map`** for NAME→ACC,DESC. |
| `scripts/filter_length1_hmms.py` | Filters a HMMER3 HMM file to remove models with `LENG 1` (incompatible with nail); appends removed model names to `data/install.log`. Used by `install.sh` for Superfamily. |
| `data/install.log` | Append-only log written by install steps (e.g. which HMM models were removed due to model_length=1). |
| `nailscan.single.sh` | Runs nail on one FASTA; by default all DBs from config; writes one merged TSV (requires **`-b`**). Options: `-h`, `--applications`, `--iprlookup`, `--goterms`. |
| `nailscan.batch.sh` | Combines FASTAs in a directory, runs nail per DB, merges per genome; writes **`<output_base>/<genome_id>.tsv.gz`** (requires **`-b`**). |
| `scripts/download_interpro_mappings.sh` | Downloads interpro2go and interpro.xml.gz, builds `data/signature2interpro.tsv` for `--iprlookup` / `--goterms`. (Also run automatically by `install.sh` if files are missing.) |
| `scripts/add_ipr_go.py` | Adds InterPro and/or GO columns to a TSV using `data/` mapping files. For PANTHER, always applied: strips `.orig.30.pir` suffix from model IDs, adds family descriptions from `data/panther/PANTHER*_HMM_classifications`, and applies best-hit filtering (one PANTHER family per protein, matching InterProScan behaviour). |
| `data/panther/PANTHER*_HMM_classifications` | PANTHER family-level name/description lookup downloaded by `install.sh` from data.pantherdb.org. Used by `add_ipr_go.py`. |

## Contributors

- **Elly Poretsky (@eporetsky)**
- **Contributions and collaborations are welcome**
