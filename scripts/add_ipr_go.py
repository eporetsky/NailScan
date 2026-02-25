#!/usr/bin/env python3
"""
Add InterPro and/or GO columns to a NailScan TSV using data/signature2interpro.tsv and data/interpro2go.
Also applies database-specific post-processing to align output with InterProScan conventions:

  PANTHER     - ID cleanup (strip .orig.30.pir suffix), best-hit filter per protein, DESC from
                data/panther/PANTHER*_HMM_classifications.
  HAMAP       - Best-hit filter per protein (highest-scoring family match only, matching
                InterProScan behaviour where each protein gets at most one HAMAP annotation).
  SUPERFAMILY - ID format: "143243.1" -> "SSF143243" (add SSF prefix, strip version suffix).
                E-value threshold: keep only hits with e-value <= 1e-4 (InterProScan default).
  Pfam        - GA domain-threshold filter (reads data/pfam/pfam_a.ga built by install.sh).
                ACC version stripped: "PF00329.26" -> "PF00329".
  NCBIfam     - TC sequence-threshold filter (reads data/ncbifam/ncbifam.tc built by install.sh).
                ACC used as canonical ID (not NAME/PRK model name); version stripped:
                "NF009141.0" -> "NF009141".

Usage: add_ipr_go.py --tsv FILE --db DBNAME --data-dir DIR [--iprlookup] [--goterms]
"""
import argparse
import csv
import glob
import os
import sys


# ── Data loaders ────────────────────────────────────────────────────────────

def load_signature2interpro(path):
    m = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            db, sig, ipr = parts[0], parts[1], parts[2]
            key = (db, sig)
            if key not in m:
                m[key] = []
            m[key].append(ipr)
    return m


def load_interpro2go(path):
    m = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("!"):
                continue
            if "InterPro:" not in line or "GO:" not in line:
                continue
            try:
                ipr_part, rest = line.split(">", 1)
                ipr = ipr_part.replace("InterPro:", "").strip().split()[0]
                go_part = rest.split(";")[-1].strip()
                if go_part.startswith("GO:"):
                    if ipr not in m:
                        m[ipr] = []
                    m[ipr].append(go_part)
            except Exception:
                continue
    return m


def load_panther_names(data_dir):
    candidates = sorted(glob.glob(os.path.join(data_dir, "panther", "PANTHER*_HMM_classifications")))
    if not candidates:
        return {}
    names = {}
    with open(candidates[-1]) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            pthr_id, name = parts[0].strip(), parts[1].strip()
            if ":" not in pthr_id and name and name != "FAMILY NOT NAMED":
                names[pthr_id] = name
    return names


def load_thresholds(path):
    """Load NAME -> (seq_thresh, dom_thresh) from a .ga or .tc file."""
    if not path or not os.path.isfile(path):
        return {}
    m = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                try:
                    m[parts[0]] = (float(parts[1]), float(parts[2]))
                except ValueError:
                    pass
    return m


# ── ID helpers ───────────────────────────────────────────────────────────────

def pthr_base_id(raw):
    """PTHR16038.orig.30.pir -> PTHR16038"""
    if not raw:
        return raw
    base = raw.split(".")[0]
    return base if base.startswith("PTHR") else raw


def ssf_id(raw_acc):
    """143243.1 -> SSF143243  (already SSF* -> unchanged)"""
    if not raw_acc:
        return raw_acc
    if raw_acc.startswith("SSF"):
        return raw_acc
    base = raw_acc.split(".")[0]
    return f"SSF{base}" if base.isdigit() else raw_acc


def strip_version(acc):
    """PF00329.26 -> PF00329 ; NF009141.0 -> NF009141"""
    if not acc:
        return acc
    dot = acc.rfind(".")
    if dot == -1:
        return acc
    suffix = acc[dot + 1:]
    # Only strip if the suffix is purely numeric (a version number)
    return acc[:dot] if suffix.isdigit() else acc


# ── Best-hit helper ──────────────────────────────────────────────────────────

def best_hit_per_target(rows):
    """Keep the row with the highest bit score for each target protein."""
    best = {}
    for row in rows:
        target = row.get("target", "")
        try:
            score = float(row.get("score") or 0)
        except (ValueError, TypeError):
            score = 0.0
        if target not in best or score > best[target][0]:
            best[target] = (score, row)
    return [row for _, (_, row) in sorted(best.items(), key=lambda x: x[0])]


# ── IPR / GO annotation ──────────────────────────────────────────────────────

def annotate_ipr_go(rows, db_xml, sig2ipr, ipr2go, do_ipr, do_go):
    for row in rows:
        acc = (row.get("ACC") or row.get("NAME") or "").strip()
        ipr_list = sig2ipr.get((db_xml, acc), [])
        if not ipr_list and acc and "." in acc:
            ipr_list = sig2ipr.get((db_xml, acc.split(".")[0]), [])
        if do_ipr:
            row["InterPro"] = ",".join(ipr_list) if ipr_list else ""
        if do_go:
            go_set = set()
            for ipr in ipr_list:
                go_set.update(ipr2go.get(ipr, []))
            row["GO"] = ",".join(sorted(go_set)) if go_set else ""


# ── DB name mapping ──────────────────────────────────────────────────────────

DB_TO_XML = {
    "pfam": "PFAM",
    "ncbifam": "NCBIFAM",
    "superfamily": "SSF",
    "hamap": "HAMAP",
    "sfld": "SFLD",
    "pirsf": "PIRSF",
    "pirsr": "PIRSR",
    "panther": "PANTHER",
    "cath": "CATHGENE3D",
    "antifam": "ANTIFAM",
}


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--iprlookup", action="store_true")
    ap.add_argument("--goterms", action="store_true")
    args = ap.parse_args()

    db = args.db.lower()
    data_dir = args.data_dir
    db_xml = DB_TO_XML.get(db, args.db.upper())

    # Databases with custom post-processing logic
    CUSTOM_DBS = {"panther", "hamap", "superfamily", "pfam", "ncbifam"}
    needs_custom = db in CUSTOM_DBS

    # Pass-through for dbs with no custom logic and no annotation requested
    if not args.iprlookup and not args.goterms and not needs_custom:
        with open(args.tsv) as f:
            sys.stdout.write(f.read())
        return

    # Load IPR/GO data if needed
    sig2ipr, ipr2go = {}, {}
    if args.iprlookup or args.goterms:
        sig2ipr_path = os.path.join(data_dir, "signature2interpro.tsv")
        ipr2go_path = os.path.join(data_dir, "interpro2go")
        if not os.path.isfile(sig2ipr_path):
            sys.stderr.write(f"Error: {sig2ipr_path} not found.\n")
            sys.exit(1)
        if args.goterms and not os.path.isfile(ipr2go_path):
            sys.stderr.write(f"Error: {ipr2go_path} not found.\n")
            sys.exit(1)
        sig2ipr = load_signature2interpro(sig2ipr_path)
        ipr2go = load_interpro2go(ipr2go_path) if args.goterms else {}

    with open(args.tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames)
        if args.iprlookup and "InterPro" not in fieldnames:
            fieldnames.append("InterPro")
        if args.goterms and "GO" not in fieldnames:
            fieldnames.append("GO")

        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")

        # ── PANTHER ──────────────────────────────────────────────────────────
        if db == "panther":
            panther_names = load_panther_names(data_dir)
            all_rows = list(reader)
            for row in all_rows:
                base_id = pthr_base_id(row.get("NAME", ""))
                row["NAME"] = base_id
                row["ACC"] = base_id
                if "DESC" in fieldnames and not (row.get("DESC") or "").strip():
                    row["DESC"] = panther_names.get(base_id, "")
            rows = best_hit_per_target(all_rows)
            annotate_ipr_go(rows, db_xml, sig2ipr, ipr2go, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(rows)
            return

        # ── HAMAP ─────────────────────────────────────────────────────────────
        # HAMAP has no GA/TC thresholds in its HMM files; InterProScan reports at
        # most one HAMAP family per protein (the best-scoring match).
        if db == "hamap":
            all_rows = list(reader)
            rows = best_hit_per_target(all_rows)
            annotate_ipr_go(rows, db_xml, sig2ipr, ipr2go, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(rows)
            return

        # ── SUPERFAMILY ───────────────────────────────────────────────────────
        # InterProScan uses E < 1e-4 for Superfamily; converts numeric ACC to SSFxxxx.
        if db == "superfamily":
            all_rows = list(reader)
            kept = []
            for row in all_rows:
                try:
                    evalue = float(row.get("evalue") or 1)
                except (ValueError, TypeError):
                    evalue = 1.0
                if evalue > 1e-4:
                    continue
                clean = ssf_id(row.get("ACC") or row.get("NAME") or "")
                row["NAME"] = clean
                row["ACC"] = clean
                kept.append(row)
            annotate_ipr_go(kept, db_xml, sig2ipr, ipr2go, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(kept)
            return

        # ── Pfam ──────────────────────────────────────────────────────────────
        # Filter by GA domain threshold; strip version from ACC.
        if db == "pfam":
            ga = load_thresholds(os.path.join(data_dir, "pfam", "pfam_a.ga"))
            all_rows = list(reader)
            kept = []
            for row in all_rows:
                name = row.get("NAME", "")
                try:
                    score = float(row.get("score") or 0)
                except (ValueError, TypeError):
                    score = 0.0
                if name in ga:
                    dom_ga = ga[name][1]
                    if score < dom_ga:
                        continue
                # Strip version from ACC
                acc = strip_version(row.get("ACC") or "")
                row["ACC"] = acc
                kept.append(row)
            annotate_ipr_go(kept, db_xml, sig2ipr, ipr2go, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(kept)
            return

        # ── NCBIfam ───────────────────────────────────────────────────────────
        # Filter by TC sequence threshold; use ACC as canonical ID (strip version).
        # The NAME field contains the internal PRK model name; the ACC field holds
        # the public NF/TIGR identifier. Use ACC for the output.
        if db == "ncbifam":
            tc = load_thresholds(os.path.join(data_dir, "ncbifam", "ncbifam.tc"))
            all_rows = list(reader)
            kept = []
            for row in all_rows:
                name = row.get("NAME", "")
                try:
                    score = float(row.get("score") or 0)
                except (ValueError, TypeError):
                    score = 0.0
                if name in tc:
                    seq_tc = tc[name][0]
                    if score < seq_tc:
                        continue
                # Promote ACC to NAME; strip version suffix
                raw_acc = row.get("ACC") or row.get("NAME") or ""
                clean_acc = strip_version(raw_acc)
                row["NAME"] = clean_acc
                row["ACC"] = clean_acc
                kept.append(row)
            annotate_ipr_go(kept, db_xml, sig2ipr, ipr2go, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(kept)
            return

        # ── All other databases ───────────────────────────────────────────────
        writer.writeheader()
        for row in reader:
            acc = (row.get("ACC") or "").strip()
            if not acc:
                acc = (row.get("NAME") or "").strip()
            ipr_list = sig2ipr.get((db_xml, acc), [])
            if not ipr_list and acc and "." in acc:
                ipr_list = sig2ipr.get((db_xml, acc.split(".")[0]), [])
            if args.iprlookup:
                row["InterPro"] = ",".join(ipr_list) if ipr_list else ""
            if args.goterms:
                go_set = set()
                for ipr in ipr_list:
                    go_set.update(ipr2go.get(ipr, []))
                row["GO"] = ",".join(sorted(go_set)) if go_set else ""
            writer.writerow(row)


if __name__ == "__main__":
    main()
