# NailScan

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

NailScan uses [nail](https://github.com/TravisWheelerLab/nail) (an Alignment Inference tooL) to run fast profile HMM search against HMM-based database such as the [Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/) database. Nail uses MMseqs2 for seeding and approximates the HMMER3 Forward/Backward algorithm for speed. NailScan is **not** meant to reproduce identical results to InterProScan or to replace it; it is intended as an alternative for HMM-based analyses where speed and scalability are priorities.

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

2. **Install nail and the Pfam database**

  ```bash
  ./install.sh
  ```

   This builds nail and, using `config.json`, downloads the Pfam database from InterPro/EBI, verifies the MD5 checksum, and unpacks it into `data/pfam/`. If `data/pfam/pfam_a.hmm` already exists, the download is skipped.

3. **Run Pfam search**

   - **Single FASTA** (requires `-f`)  
     `./nailscan.single.sh -f path/to/sequences.fasta`  
     Writes `results/<basename>.pfam.tsv` (tab-separated, with optional Pfam ACC/DESC columns).

   - **Batch (directory of FASTA files)**  
     `./nailscan.batch.sh -f path/to/fasta_dir`  
     Collects all `.fa`, `.faa`, and `.fasta` files in the directory, runs one nail search, then splits results into `results/pfam/<genome_id>.tsv.gz` (default output).

### Optional: test run

The `test/` directory contains three small FASTA files (`test1.fa`, `test2.fa`, `test3.fa`) with a few dozen sequences each. After installation, you can try:

**Single FASTA** (one file):

```bash
./nailscan.single.sh -f test/test1.fa -t 72
# → writes results/test1.pfam.tsv
```

**Batch** (all FASTA files in `test/`; output defaults to results/pfam/):

```bash
./nailscan.batch.sh -f test/ -t 72
# → writes results/pfam/test1.tsv.gz, results/pfam/test2.tsv.gz, results/pfam/test3.tsv.gz
```

## Threading and resources

- Use **`-t N`** to set the number of threads (default: all cores via `nproc`, or 16).

  ```bash
  ./nailscan.single.sh -f my.faa -t 72
  ./nailscan.batch.sh -f /path/to/fastas -t 72
  ```

- **Memory:** Batch mode can use a lot of RAM. As a rough guide: a run with **~5.5 million sequences** used **~240 GB RAM** and completed in **~25 minutes** with **72 threads**. Plan resources accordingly for large batches.

## Default paths

| Script              | Pfam HMM (`-h`)      | Input (`-f`)     | Output (`-o`)      |
|---------------------|----------------------|------------------|---------------------|
| `nailscan.single.sh`| `data/pfam/pfam_a.hmm` | **required**    | `results/`          |
| `nailscan.batch.sh` | `data/pfam/pfam_a.hmm` | `fasta/`        | `results/pfam/`     |

Options: **`-f`** fasta file (single) or directory (batch), **`-o`** output file/dir, **`-h`** Pfam HMM path, **`-t`** threads. Use `--help` for usage.

```bash
./nailscan.single.sh -f <input.fasta> [-h pfam_hmm] [-o output_dir] [-t N]
./nailscan.batch.sh [-f fasta_dir] [-h pfam_hmm] [-o results_dir] [-t N]
```

## Output format

- **nailscan.single.sh** writes e.g. `results/<basename>.pfam.tsv`: tab-separated, one header row. Columns: `target`, `NAME`, `ACC`, `DESC`, `target_start`, `target_end`, `query_start`, `query_end`, `score`, `bias`, `evalue`, `cell_frac`.

- **nailscan.batch.sh** writes `results/pfam/<genome_id>.tsv.gz`: same column layout, one gzipped TSV per input file/genome.

**E‑value filtering:** Pfam/nail hits may include weak matches. Consider filtering by E‑value (and/or score) in downstream analysis.

## Configuration

- **`config.json`** — Used by `install.sh` to download the Pfam database. Example:

  ```json
  {
    "pfam": {
      "url": "https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/6/6.0/pfam/pfam-38.1.tar.gz",
      "md5_url": "https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/6/6.0/pfam/pfam-38.1.tar.gz.md5"
    }
  }
  ```

  The installer checks the MD5 before unpacking. If `data/pfam/pfam_a.hmm` is already present, the download step is skipped.

## What each file does

| File                 | Purpose |
|----------------------|--------|
| `environment.yml`     | Conda env with **Rust** (to build nail) and **MMseqs2** (required by nail at runtime). |
| `config.json`        | Pfam download URL and MD5 URL for `install.sh`. |
| `install.sh`         | Clones nail, builds with Cargo, installs `nail` into the env; optionally downloads Pfam to `data/pfam/` using `config.json` (with MD5 check). Defines **`make_pfam_map`** for NAME→ACC,DESC mapping. |
| `nailscan.single.sh` | Runs nail on one FASTA; writes one TSV (with optional ACC/DESC from Pfam map). |
| `nailscan.batch.sh`  | Combines all `.fa`/`.faa`/`.fasta` in a directory, runs one nail search, splits results to gzipped TSVs in `results/pfam/`. |
| `scripts/make_pfam_map.sh` | Parses a Pfam HMM and writes NAME, ACC, DESC map; used by the scan scripts when the map is missing. |

## Contributors
- **Elly Poretsky (@eporetsky)**
- **Contributions and collaborations are welcome**
