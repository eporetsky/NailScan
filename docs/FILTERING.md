# NailScan — Database Filtering and Output Format Reference

This document describes the post-processing, filtering, and output formatting
decisions applied by NailScan to each supported database.  It records the
rationale behind every design choice so it is easy to revisit or adjust
individual databases in `scripts/add_ipr_go.py` and `config.json`.

---

## 1. Output columns

| Column | Description |
|--------|-------------|
| `target` | Protein / gene identifier (query sequence) |
| `Analysis` | Database name (canonical capitalisation — see §3) |
| `ACC` | Canonical accession for the matched model |
| `DESC` | Human-readable description of the model |
| `target_start` | Start position of the match on the protein (1-based) |
| `target_end` | End position of the match on the protein (inclusive) |
| `query_start` | Start of the alignment on the HMM |
| `query_end` | End of the alignment on the HMM |
| `score` | Bit score from nail |
| `bias` | Composition bias |
| `evalue` | E-value |
| `cell_frac` | Fraction of DP matrix cells computed by nail |
| `InterPro` | InterPro entry ID (requires `--iprlookup`) |
| `IPR_desc` | Short description of the InterPro entry (requires `--iprlookup`) |
| `GO` | GO term IDs, semicolon-separated (requires `--goterms`) |

The `NAME` column (raw HMM model name) is present in intermediate per-database
temporary files used internally by `add_ipr_go.py` but is **stripped from the
final merged output**.  This was a deliberate design choice: downstream tools
should use `ACC` as the canonical identifier.

Columns are sorted by `target` (gene ID) in the final TSV.

---

## 2. Database status

| Database | Enabled | Filtering | Notes |
|----------|---------|-----------|-------|
| **Pfam** | yes | GA domain-level threshold | §4.1 |
| **Gene3D** | yes | Best-hit per G3DSA domain; e-value ≤ 1e-4 | §4.2 |
| **SUPERFAMILY** | yes | E-value ≤ 1e-4; non-overlapping hit selection | §4.3 |
| **PANTHER** | yes | Best-hit per protein | §4.4 |
| **NCBIfam** | yes | TC sequence-level threshold | §4.5 |
| SFLD | **disabled** | — | §5.4 |
| **AntiFam** | yes | No extra filtering (flags spurious ORFs) | §4.7 |
| Hamap | **disabled** | — | §5.1 |
| PIRSF | **disabled** | — | §5.2 |
| PIRSR | **disabled** | — | §5.3 |

Databases not supported at all (format incompatibility):
Prosite, SMART, CDD, PRINTS, Phobius, TMHMM.
See `docs/INTERPRO_DATABASES.md` for details.

---

## 3. Canonical database labels

`nailscan.single.sh` maps raw database directory names to the canonical
capitalisation used in the `Analysis` column:

| Raw name | `Analysis` label |
|----------|-----------------|
| `pfam` | `Pfam` |
| `ncbifam` | `NCBIfam` |
| `cath` | `Gene3D` |
| `superfamily` | `SUPERFAMILY` |
| `panther` | `PANTHER` |
| `sfld` | `SFLD` |
| `antifam` | `AntiFam` |
| `hamap` | `Hamap` |
| `pirsf` | `PIRSF` |
| `pirsr` | `PIRSR` |

---

## 4. Enabled databases — filtering details

### 4.1 Pfam

**Source files:** `data/pfam/pfam_a.hmm`, `data/pfam/pfam_a.ga`,
`data/pfam/pfam_a.dat`

**Filtering (three steps):**

1. **Domain-level GA threshold** — each hit must have `score ≥ GA_dom` for
   its model.  Cutoffs extracted into `pfam_a.ga` at install time.

2. **Sequence-level GA threshold** — for each `(protein, model)` pair, all
   domain scores are summed as a proxy for HMMER's sequence-level score.  If
   the sum is below `GA_seq`, all domains for that pair are dropped.  For most
   Pfam models `GA_seq == GA_dom`, so this only affects models where a single
   weak domain passes the domain cutoff but the protein lacks cumulative
   evidence.

3. **Clan-based overlap deduplication** — `pfam_a.dat` contains `#=GF CL`
   clan assignments grouping related Pfam families (e.g. all kinase subtypes
   share a clan).  Within each protein, if two hits from the **same clan**
   overlap by more than 50 % of the shorter hit's length, the lower-scoring
   one is removed.  The 50 % threshold (vs. a more aggressive 30 %) avoids
   discarding genuinely co-occurring adjacent same-clan domains (e.g. tandem
   repeats, multi-domain architectures) which caused false negatives.

- ACC version suffix stripped: `PF00329.26` → `PF00329`.

---

### 4.2 Gene3D (CATH)

**Source files:** `data/cath/gene3d_main.hmm`,
`data/cath/model_to_family_map.tsv`

**Filtering:**
1. E-value ≤ 1e-4.
2. Bitscore ≥ 10 (matching `cath-resolve-hits --worst-permissible-bitscore 10`).
3. Map each HMM model name to its `G3DSA:x.x.x.x` domain via
   `model_to_family_map.tsv` (model name split on `-`, first token looked up).
4. **Global per-protein non-overlapping selection** (overlap threshold 20%):
   all Gene3D hits for a protein are sorted by score and the greedy
   non-overlapping algorithm is applied across all G3DSA domains together.

**Rationale:** InterProScan uses `cath-resolve-hits`, a binary that solves
the domain architecture problem as a dynamic-programming optimisation:
maximise the total bit score of a non-overlapping set of domain hits.
Default filters: e-value ≤ 0.001, bitscore ≥ 10, 10-residue boundary trim.
The previous NailScan approach (`best_hit_per_domain`) only removed redundancy
within the same G3DSA domain but left cross-domain overlapping hits intact,
causing FPs when two different G3DSA families matched the same protein region.
Switching to global per-protein non-overlapping selection resolves these
cross-domain conflicts; the greedy approximation is stricter than the DP
solver at boundaries (20% threshold vs ~10 residue trim) to compensate.

**Note on `cath-resolve-hits`:** The binary ships with InterProScan and is
available at `interproscan-5.77-108.0-64-bit/bin/gene3d/4.3.0/cath-resolve-hits`.
A future improvement could convert nail's tbl-out to the `raw_with_scores`
format and pipe it through `cath-resolve-hits` directly for fully faithful
DP-based domain architecture resolution.

---

### 4.3 SUPERFAMILY

**Source files:** `data/superfamily/hmmlib_1.75`
(converted to HMMER3 p7 ASCII during install; COMPO lines inserted)

**Filtering:**
1. E-value ≤ **1e-3** (see rationale below).
2. ID format conversion: `143243.1` → `SSF143243`.
3. **Global per-protein non-overlapping hit selection** (overlap threshold
   35%, matching `$percentsame = 0.35` in InterProScan's `ass3.pl`):
   - All SSF hits for a protein are sorted by score (descending) together,
     regardless of which SSF family they belong to.
   - A hit is rejected if it overlaps an already-accepted hit (from **any**
     SSF family) by more than 35 % of the shorter hit's length.

**Rationale (e-value):** `ass3.pl` states `0.0001` as the default cutoff, but
computes an adjusted `tempcutoff = cutoff / (n_sfs / totmods)` where `n_sfs`
is the number of unique superfamilies matched by the input protein and `totmods`
is the total number of HMM models.  For a typical protein matching a few dozen
superfamilies out of ~2,000+ models, this ratio is <<1, making the effective
cutoff 10–100× more permissive than 1e-4.  Using 1e-3 approximates this
adjustment and improves recall without sacrificing precision (the global
non-overlapping filter handles precision).

**Rationale (overlap filter):** The filter is applied globally, not per SSF
family.  InterProScan's `ass3.pl` processes all superfamily hits together: a
hit from SSF456 that overlaps a higher-scoring SSF123 hit by > 35 % is still
rejected, even though they are different superfamilies.  An earlier NailScan
implementation grouped hits by `(protein, domain_acc)` before filtering, which
left all cross-family overlaps uncollapsed and was the main cause of SUPERFAMILY
over-prediction.

**HMM format note:** The SUPERFAMILY HMM file ships in HMMER3 format but
lacks `COMPO` lines that nail's parser requires.  The install script uses
`awk` to insert synthetic `COMPO` lines (equal amino-acid frequencies) so
nail can parse the file without any source-level changes.

---

### 4.4 PANTHER

**Source files:** `data/panther/famhmm/binHmm` (directory of per-family HMMs),
`data/panther/PANTHER*_HMM_classifications`

**Filtering:**
1. ID cleanup: strip suffix from HMM model filenames
   (e.g. `PTHR22842.orig.30.pir` → `PTHR22842`).
2. **Best-hit per protein:** keep only the single top-scoring PANTHER family
   for each protein.
3. DESC populated from `PANTHER*_HMM_classifications` using the cleaned
   family ID.

**Rationale:** InterProScan reports one PANTHER family per protein (the
best-matching family in the tree).  Raw HMM search returns the 50–75 closest
families across the entire PANTHER tree, most of which are false positives for
a best-classification task.

---

### 4.5 NCBIfam

**Source files:** `data/ncbifam/ncbifam.hmm`, `data/ncbifam/ncbifam.tc`

**Filtering:**
- Apply the **TC (Trusted Cutoff)** sequence-level threshold from each NCBIfam
  HMM, extracted into `ncbifam.tc` at install time.
- Only hits with `score ≥ TC_sequence` for that model are kept.
- ACC version suffix stripped: `NF009141.0` → `NF009141`.

**Rationale:** NCBIfam curators define TC thresholds to distinguish true
members from noise.  Without them, nail returns many borderline hits.
InterProScan applies the same TC cutoffs.

---

### ~~4.6 SFLD~~ — moved to §5.4 (disabled)

---

### 4.7 AntiFam

**Source files:** `data/antifam/AntiFam.hmm`

**Filtering:** None beyond nail's own e-value output.

AntiFam models flag spurious ORFs and translation artefacts.  Hits to any
AntiFam model are significant regardless of score ranking.

---

## 5. Disabled databases

### 5.1 Hamap

`"enabled": false` in `config.json`.

**Reason:** Hamap uses **Prosite-matrix normalised scores (N_SCORE)**, not
HMM bit scores, for per-profile thresholding.  The cutoffs are stored in
`hamap.prf` and have no direct equivalent in HMM bit-score space.  When
searched with nail using HMM bit scores, Hamap produces approximately 10×
more predictions than InterProScan (benchmark precision ≈ 0.10 vs. expected
≈ 1.0 for a well-tuned database).

**To re-enable:** `--appl hamap` (pass-through to both single and batch
scripts).  Filtering improvements would require either: (a) deriving
per-model bit-score cutoffs empirically from a reference proteome, or (b)
running Hamap through its native PfsearchV3 engine and merging the results.

---

### 5.2 PIRSF

`"enabled": false` in `config.json`.

**Reason:** PIRSF families are **hierarchical**: a broad parent family (e.g.
"Protein kinase") may have hundreds of more specific child families.
InterProScan reports only the most specific matching child family using
length, coverage, and score criteria defined in `pirsf.dat`.  Without the
full hierarchy-aware scoring algorithm, plain HMM search reports both parent
and child families for most proteins, yielding ~4× the number of true
positives (benchmark precision ≈ 0.25).

A partial hierarchy filter is implemented in `add_ipr_go.py`
(`filter_pirsf_hierarchy`) but it does not fully replicate InterProScan's
length/coverage checks, leading to residual over-prediction.

**To re-enable:** `--appl pirsf`.

---

### 5.4 SFLD

`"enabled": false` in `config.json`.

**Reason:** SFLD requires two levels of filtering that cannot both be applied
with nail alone:

1. **GA domain-level threshold** — implemented and working.

2. **Per-residue active-site matching** — InterProScan's `sfld_postprocess`
   binary reads HMMER's Stockholm alignment output (`hmmsearch -A`) and
   cross-checks each hit against the catalytic residue patterns in
   `sfld_sites.annot`.  Proteins that pass GA but lack the required catalytic
   residues are discarded.  Nail does not produce HMMER-compatible alignment
   files, so this filter cannot be applied; it causes ~20–25% over-prediction.

3. **Hierarchy filter (SFLDS → SFLDG → SFLDF)** — an attempt to apply the
   same "keep most specific" rule as PIRSF was made using `sfld_hierarchy.tsv`,
   but this caused recall to drop from ~0.82 to ~0.41.  The SFLD hierarchy
   file contains redundant rows and the TSV does not cleanly represent a
   one-to-one parent–child tree, so the filter incorrectly removed many
   legitimate group-level (SFLDG) hits that have no matching family-level
   (SFLDF) entry.  The hierarchy filter code is preserved in
   `add_ipr_go.py` (`load_sfld_hierarchy`, `filter_sfld_hierarchy`) but is
   no longer called.

**To re-enable:** `--appl sfld`.  Results will use only the GA threshold; the
residue-level filter requires HMMER's full output files.

**Future path:** Run `hmmsearch --cut_ga -A sfld.aln --domtblout sfld.dom`
for the SFLD HMM database, then pipe through `scripts/sfld_postprocess` (the
`sfld_postprocess` binary from InterProScan is included in `scripts/`).  This
would provide fully accurate results at the cost of a separate, slower HMMER
search for this database only.

### 5.3 PIRSR

`"enabled": false` in `config.json`.

**Reason:** PIR Site Rules annotate **specific functional residues** (e.g.
catalytic sites, metal-binding sites).  They require positional evidence —
matching which residues within an alignment are conserved — not just a global
HMM score.  A plain `nail search` produces hit boundaries, not residue-level
annotations, so PIRSR predictions are conceptually unsound and cannot be
validated against InterProScan output.

**To re-enable:** `--appl pirsr`.

---

## 6. InterPro / GO annotation (--iprlookup, --goterms)

Required data files (generated or downloaded by `install.sh` and
`scripts/download_interpro_mappings.sh`):

| File | Description |
|------|-------------|
| `data/signature2interpro.tsv` | member-DB accession → InterPro ID |
| `data/interpro2go` | InterPro ID → GO terms |
| `data/interpro_names.tsv` | InterPro ID → short description |

`data/interpro_names.tsv` is built during `install.sh` by parsing
`data/interpro.xml.gz`.  If `interpro.xml.gz` is absent, the
`IPR_desc` column is omitted (no error).

---

## 7. Benchmark results (AtCol-0 Arabidopsis thaliana Col-0, 35,386 proteins)

Benchmark run against InterProScan 5.77 output using `scripts/benchmark.py`.
Enabled databases only; Hamap and PIRSF excluded.

| Database | Precision | Recall | F1 |
|----------|-----------|--------|----|
| Pfam | ~0.94 | ~0.88 | ~0.91 |
| Gene3D | ~0.88 | ~0.72 | ~0.79 |
| SUPERFAMILY | ~0.91 | ~0.85 | ~0.88 |
| PANTHER | ~0.97 | ~0.91 | ~0.94 |
| NCBIfam | ~0.96 | ~0.89 | ~0.92 |
| SFLD | ~0.90 | ~0.84 | ~0.87 |

These figures reflect the benchmark used to guide filtering design decisions.
Re-run `scripts/benchmark.py` against updated data to get current numbers.

---

## 8. Architecture notes

### Why intermediate TSVs still contain NAME

`nailscan.single.sh` (and formerly `nailscan.batch.sh`) write per-database
temporary TSVs that include a `NAME` column (the raw HMM model name from nail's
`--tbl-out`).  `add_ipr_go.py` relies on `NAME` for:

- **Pfam / NCBIfam / SFLD:** looking up GA/TC thresholds (the threshold files
  are keyed by model name, not accession).
- **PANTHER:** stripping the `.orig.30.pir` suffix.
- **Hamap:** using `NAME` as the fallback for `ACC` when `ACC` is empty.

The `NAME` column is stripped from the **final merged output** by
`add_ipr_go.py`.  Removing it from the intermediate files broke threshold
lookups and PANTHER ID parsing (regression observed Dec 2025).

### Batch script design

`nailscan.batch.sh` is a thin wrapper around `nailscan.single.sh`.  All HMM
search, filtering, and column manipulation logic lives exclusively in
`nailscan.single.sh` and `scripts/add_ipr_go.py`.  The batch script only:

1. Chunks a directory of genome FASTAs into batches.
2. Prefixes every protein ID with `genome_id|` so a single nail run covers
   multiple genomes without ID collisions.
3. Calls `nailscan.single.sh` on each combined batch FASTA.
4. Splits the combined output TSV back into per-genome files and compresses
   them with gzip.

This means bug fixes and new features in `nailscan.single.sh` automatically
apply to batch mode without any separate edits.

### SUPERFAMILY HMM format fix

SUPERFAMILY HMMs ship in HMMER3 text format but lack `COMPO` lines.  nail's
HMM parser requires a `COMPO` line immediately after the model header.
`install.sh` inserts synthetic `COMPO` lines (uniform composition) via `awk`
into a processed copy (`hmmlib_1.75.nail`) during installation.  This avoids
any modification to nail's source code.  The processed file is what nail
searches; the original is preserved as-is.

### FunFam

Gene3D/CATH releases a FunFam sub-classification alongside the main structural
domain models.  InterProScan reports FunFam assignments in addition to Gene3D
domain hits.  NailScan does not currently run FunFam because the `Gene3D` label
in the `Analysis` column already refers to structural domain-level hits.
FunFam would require a separate HMM library and its own filtering.

---

## 9. Re-enabling or adding a database

1. Set `"enabled": true` in `config.json` (or omit the key; default is true).
2. Add a filtering block to `scripts/add_ipr_go.py` in the `CUSTOM_DBS`
   dispatch dict.
3. Add the database to `_SKIP_DBS` in `scripts/benchmark.py` until filtering
   is validated, then remove it.
4. Run `./install.sh` to download the HMM library.
5. Benchmark with `scripts/benchmark.py` against an InterProScan reference.
