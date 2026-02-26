#!/usr/bin/env python3
"""
Compare NailScan output against InterProScan (or another NailScan run) at the
(protein, database, accession) level.

Usage:
  python3 scripts/benchmark.py --nailscan results/myrun.tsv --ipr test/AtCol-0.fa.tsv
  python3 scripts/benchmark.py --nailscan results/run_a.tsv --nailscan2 results/run_b.tsv

NailScan TSV format (with Analysis column):
  Analysis  target  ACC  DESC  target_start ...  evalue ...

InterProScan TSV format (no header):
  protein_id  md5  len  db  acc  desc  start  end  score  T/F  date  ipr_id  ipr_desc  go  pathway

Output:
  Per-database precision / recall / F1 table, plus overall summary.
  Optionally write per-protein detail to a TSV with --detail out.tsv.
"""
import argparse
import csv
import sys
from collections import defaultdict


# ── ID normalisation ─────────────────────────────────────────────────────────

def strip_version(acc: str) -> str:
    """PF00329.26 -> PF00329 ; NF009141.0 -> NF009141"""
    if not acc:
        return acc
    dot = acc.rfind(".")
    if dot == -1:
        return acc
    suffix = acc[dot + 1:]
    return acc[:dot] if suffix.isdigit() else acc


# Map common DB name variants to a canonical form
_DB_ALIASES = {
    "pfam": "Pfam", "PFAM": "Pfam",
    "ncbifam": "NCBIfam", "NCBIFAM": "NCBIfam",
    "superfamily": "SUPERFAMILY",
    "hamap": "Hamap", "HAMAP": "Hamap",
    "panther": "PANTHER",
    "pirsf": "PIRSF",
    "pirsr": "PIRSR",
    "cath": "Gene3D", "cathgene3d": "Gene3D",
    "antifam": "AntiFam",
    "sfld": "SFLD",
    "gene3d": "Gene3D", "Gene3D": "Gene3D",
    "SUPERFAMILY": "SUPERFAMILY",
    "Hamap": "Hamap",
    "PANTHER": "PANTHER",
    "NCBIfam": "NCBIfam",
    "Pfam": "Pfam",
}

# Databases not currently run by NailScan (or whose results are not yet reliable
# enough to include in the benchmark). Commented out here so they are easy to
# re-enable individually when support is added.
_SKIP_DBS = {
    # Not in NailScan config:
    "CDD",
    "PRINTS",
    "SMART",
    "ProSitePatterns",
    "ProSiteProfiles",
    # Gene3D sub-classification, not a separate NailScan analysis:
    "FunFam",
    # Excluded from default benchmark (comment out to re-enable):
    "AntiFam",      # rarely annotates proteins; not a priority
    "PIRSR",        # PIR Site Rules require residue-level detection algorithm
    "Hamap",        # Prosite-matrix N_SCORE thresholds not transferable to HMM bit scores
    "PIRSF",        # hierarchy-aware scoring not fully implemented; over-predicts parent families
    "SFLD",         # requires per-residue active-site matching (sfld_postprocess); nail does not
                    # produce HMMER alignment files needed for it; hierarchy filter caused recall ~0.4
}


def norm_db(db: str) -> str:
    return _DB_ALIASES.get(db, db)


def norm_acc(acc: str, db: str) -> str:
    """Normalise accession: strip version, apply per-db rules."""
    acc = strip_version(acc.strip()) if acc else ""
    # Superfamily: bare number -> SSFnumber
    if db == "SUPERFAMILY" and acc and not acc.startswith("SSF"):
        base = acc.split(".")[0]
        if base.isdigit():
            acc = f"SSF{base}"
    return acc


# ── Loaders ──────────────────────────────────────────────────────────────────

def load_nailscan(path: str) -> dict:
    """Return dict (protein, db, acc) -> set of (start, end) tuples."""
    hits = defaultdict(set)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        # Accept both with and without the Analysis column, and both old
        # (NAME present) and new (NAME removed) output formats.
        has_analysis = "Analysis" in (reader.fieldnames or [])
        for row in reader:
            db_raw = row.get("Analysis", "unknown") if has_analysis else "unknown"
            db = norm_db(db_raw)
            if db in _SKIP_DBS:
                continue
            protein = row.get("target", "")
            acc_raw = row.get("ACC") or row.get("NAME") or ""
            acc = norm_acc(acc_raw, db)
            if not protein or not acc:
                continue
            start = row.get("target_start", "")
            end = row.get("target_end", "")
            hits[(protein, db, acc)].add((start, end))
    return hits


def load_ipr(path: str) -> dict:
    """
    Load InterProScan TSV (no header).
    Columns: protein_id, md5, len, db, acc, desc, start, end, score, T/F,
             date, ipr_id, ipr_desc, go, pathway
    """
    hits = defaultdict(set)
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                continue
            protein = cols[0].strip()
            db = norm_db(cols[3].strip())
            acc_raw = cols[4].strip()
            acc = norm_acc(acc_raw, db)
            start = cols[6] if len(cols) > 6 else ""
            end = cols[7] if len(cols) > 7 else ""
            if not protein or not acc or db in {"MobiDBLite", "Coils"} | _SKIP_DBS:
                continue
            hits[(protein, db, acc)].add((start, end))
    return hits


# ── Metrics ──────────────────────────────────────────────────────────────────

def compute_stats(pred: dict, ref: dict):
    """
    Compare prediction vs reference at the (protein, db, acc) key level.
    Returns dict: db -> {tp, fp, fn, precision, recall, f1}.
    """
    pred_keys = set(pred.keys())
    ref_keys = set(ref.keys())

    # Collect all databases
    all_dbs = {k[1] for k in pred_keys | ref_keys}
    stats = {}

    for db in sorted(all_dbs):
        pred_db = {k for k in pred_keys if k[1] == db}
        ref_db = {k for k in ref_keys if k[1] == db}

        tp = len(pred_db & ref_db)
        fp = len(pred_db - ref_db)
        fn = len(ref_db - pred_db)

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = (2 * precision * recall / (precision + recall)
              if (precision + recall) > 0 else 0.0)

        stats[db] = {
            "tp": tp, "fp": fp, "fn": fn,
            "pred_total": tp + fp,
            "ref_total": tp + fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        }
    return stats


def overall_stats(stats: dict):
    tp = sum(v["tp"] for v in stats.values())
    fp = sum(v["fp"] for v in stats.values())
    fn = sum(v["fn"] for v in stats.values())
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
    return {"tp": tp, "fp": fp, "fn": fn,
            "pred_total": tp + fp, "ref_total": tp + fn,
            "precision": precision, "recall": recall, "f1": f1}


# ── Detail output ─────────────────────────────────────────────────────────────

def write_detail(pred: dict, ref: dict, path: str):
    pred_keys = set(pred.keys())
    ref_keys = set(ref.keys())
    all_keys = sorted(pred_keys | ref_keys)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["protein", "db", "acc", "in_pred", "in_ref", "status"])
        for key in all_keys:
            protein, db, acc = key
            in_pred = key in pred_keys
            in_ref = key in ref_keys
            if in_pred and in_ref:
                status = "TP"
            elif in_pred:
                status = "FP"
            else:
                status = "FN"
            w.writerow([protein, db, acc, int(in_pred), int(in_ref), status])
    print(f"Detail written to {path}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Benchmark NailScan vs InterProScan or another run")
    ap.add_argument("--nailscan", required=True,
                    help="NailScan TSV (with Analysis column)")
    ap.add_argument("--ipr", default=None,
                    help="InterProScan TSV (no-header, standard format)")
    ap.add_argument("--nailscan2", default=None,
                    help="Second NailScan TSV (compare two runs)")
    ap.add_argument("--detail", default=None,
                    help="Write per-(protein,db,acc) TP/FP/FN detail to this TSV")
    ap.add_argument("--dbs", default=None,
                    help="Comma-separated list of databases to include (default: all)")
    ap.add_argument("--no-filter", action="store_true",
                    help="Do not filter reference to predicted proteins (compare vs full genome)")
    args = ap.parse_args()

    if not args.ipr and not args.nailscan2:
        ap.error("Provide --ipr or --nailscan2 as the reference")

    print(f"Loading prediction: {args.nailscan}")
    pred = load_nailscan(args.nailscan)

    if args.nailscan2:
        print(f"Loading reference (nailscan): {args.nailscan2}")
        ref = load_nailscan(args.nailscan2)
    else:
        print(f"Loading reference (InterProScan): {args.ipr}")
        ref = load_ipr(args.ipr)

    # By default filter reference to proteins that appear in the prediction,
    # so recall is meaningful for the queried subset (not the whole genome).
    if not args.no_filter:
        pred_proteins = {k[0] for k in pred}
        ref = {k: v for k, v in ref.items() if k[0] in pred_proteins}
        print(f"(Reference filtered to {len(pred_proteins)} proteins present in prediction)")

    # Filter databases if requested
    if args.dbs:
        keep = {norm_db(d.strip()) for d in args.dbs.split(",")}
        pred = {k: v for k, v in pred.items() if k[1] in keep}
        ref = {k: v for k, v in ref.items() if k[1] in keep}

    stats = compute_stats(pred, ref)
    # Exclude skipped databases from overall stats
    active_stats = {db: s for db, s in stats.items() if db not in _SKIP_DBS}
    ov = overall_stats(active_stats)

    # Print table
    col_w = 18
    hdr = (f"{'Database':<{col_w}}  {'Pred':>6}  {'Ref':>6}  "
           f"{'TP':>6}  {'FP':>6}  {'FN':>6}  "
           f"{'Prec':>7}  {'Recall':>7}  {'F1':>7}")
    sep = "-" * len(hdr)
    print()
    print(hdr)
    print(sep)
    for db, s in sorted(active_stats.items()):
        print(f"{db:<{col_w}}  {s['pred_total']:>6}  {s['ref_total']:>6}  "
              f"{s['tp']:>6}  {s['fp']:>6}  {s['fn']:>6}  "
              f"{s['precision']:>7.3f}  {s['recall']:>7.3f}  {s['f1']:>7.3f}")
    print(sep)
    print(f"{'OVERALL':<{col_w}}  {ov['pred_total']:>6}  {ov['ref_total']:>6}  "
          f"{ov['tp']:>6}  {ov['fp']:>6}  {ov['fn']:>6}  "
          f"{ov['precision']:>7.3f}  {ov['recall']:>7.3f}  {ov['f1']:>7.3f}")

    if args.detail:
        write_detail(pred, ref, args.detail)


if __name__ == "__main__":
    main()
