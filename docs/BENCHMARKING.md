# NailScan benchmarking

This document describes how `scripts/benchmark.py` compares NailScan output to InterProScan (or a second NailScan run). It clarifies what is and is not measured so results can be interpreted correctly and extended later.

---

## 1. What the benchmark compares

The script evaluates agreement at the **(protein, database, accession)** level:

- **Key:** `(protein_id, db_name, domain_accession)` — e.g. `(AT1G01010, Pfam, PF00001)`.
- **Metric:** For each such key, did both the prediction and the reference report that domain on that protein?

So the question answered is: *“How often does NailScan call the same domain types per protein as the reference?”*

Precision, recall, and F1 are computed per database and overall from counts of:

- **TP** — key present in both prediction and reference  
- **FP** — key present only in prediction  
- **FN** — key present only in reference  

---

## 2. What the benchmark does *not* compare

- **Locations / overlap:** Start and end coordinates are read from both TSVs and stored, but **never used** in the comparison. If NailScan predicts Pfam PF00001 on protein X at 10–100 and InterProScan has PF00001 on X at 500–600, it still counts as one TP. The benchmark does not check whether predicted regions overlap.
- **Number of instances:** Multiple hits of the same domain on the same protein (e.g. two copies of PF00001) are collapsed into a single key. The script does not compare how many times that domain appears or whether instance counts match.

So the benchmark is **presence/absence of domain annotations per protein**, not coordinate-level or instance-level agreement.

---

## 3. How multiple similar domains are handled

- Data structure: each key maps to a **set of (start, end)** tuples (one per hit).
- Statistics use only the **keys**. So for a given (protein, db, acc):
  - One hit in prediction and two in reference → one TP (one key match).
  - Two hits in prediction and one in reference → one TP.
- Multi-copy or repeated domains on the same protein are therefore treated as a single “present/absent” decision; the benchmark does not reward or penalise agreement on the number of copies.

---

## 4. Running the benchmark

```bash
# NailScan vs InterProScan
python3 scripts/benchmark.py --nailscan results/myrun.tsv --ipr results/reference.ipr.tsv

# Restrict to specific databases
python3 scripts/benchmark.py --nailscan results/myrun.tsv --ipr results/reference.ipr.tsv --dbs pfam,superfamily

# NailScan vs NailScan (e.g. two runs or parameter sets)
python3 scripts/benchmark.py --nailscan results/run_a.tsv --nailscan2 results/run_b.tsv

# Per-(protein, db, acc) detail TSV
python3 scripts/benchmark.py --nailscan results/myrun.tsv --ipr results/reference.ipr.tsv --detail detail.tsv
```

By default, the reference is restricted to proteins that appear in the prediction so recall is over the same set. Use `--no-filter` to compare against the full reference (e.g. full genome).

---

## 5. Validity and interpretation

- **Valid for:** Measuring how well NailScan agrees with the reference on **which domain types** are present **per protein**. That is a useful, conservative annotation-level benchmark.
- **Not valid for:** Inferring agreement on domain **positions**, **overlap**, or **number of instances** per protein.

When citing benchmark results, it is accurate to say that precision/recall/F1 reflect **domain-type presence/absence per protein**, not location or copy-number agreement.

---

## 6. Future work (to revisit later)

Possible extensions when needed:

- **Coordinate-level metrics:** Define TP as overlapping (or overlapping above a threshold) predicted vs reference regions for the same (protein, db, acc); compute instance-level precision/recall and optionally overlap-based scores (e.g. Jaccard, overlap fraction).
- **Instance count:** Optionally report or penalise agreement on the number of hits per (protein, db, acc).
- **Documentation in code:** Add a short note in the script docstring and above `compute_stats()` stating that the current implementation is presence/absence only and that coordinates are loaded but not used in the comparison.

Until then, this document is the reference for how the benchmark behaves.
