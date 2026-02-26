#!/usr/bin/env python3
"""
Post-process NailScan TSV: apply database-specific ID cleanup, score filtering,
and optionally add InterPro / GO / IPR_desc columns.

Database-specific logic
-----------------------
  PANTHER     - ID cleanup (strip .orig.30.pir suffix), best-hit per protein,
                DESC from data/panther/PANTHER*_HMM_classifications.
  HAMAP       - Best-hit per protein; ACC set from NAME when ACC is empty.
  SUPERFAMILY - ID format: "143243.1" -> "SSF143243".
                E-value threshold: <= 1e-4.
                Non-overlapping hit selection per (protein, domain):
                  hits are sorted by score and greedily accepted when their
                  overlap with already-accepted hits is < 25% of the shorter
                  hit (matching InterProScan's domain-collapsing behaviour).
  Pfam        - GA domain-threshold filter (reads data/pfam/pfam_a.ga).
                ACC version stripped: "PF00329.26" -> "PF00329".
  NCBIfam     - TC sequence-threshold filter (reads data/ncbifam/ncbifam.tc).
                ACC used as canonical ID; version stripped.
  SFLD        - GA domain-threshold filter (reads data/sfld/sfld.ga).
                Hierarchy filter: SFLD has three levels (SFLDS > SFLDG > SFLDF);
                when a protein matches both a parent and a child, only the most
                specific (child) entry is kept (data/sfld/sfld_hierarchy.tsv).
  PIRSF       - Score threshold from data/pirsf/pirsf.dat (avg - 3.5*std).
                Hierarchy filter: when a protein matches both a parent and a
                child PIRSF family, only the child (more specific) is kept.
                DESC set from HMM family name when DESC is empty.
  CATH        - Model name mapped to G3DSA:x.x.x.x via
                data/cath/model_to_family_map.tsv.
                Bitscore >= 10; e-value <= 1e-4.
                Global per-protein non-overlapping selection (overlap
                threshold 0.20) replaces the previous best-hit-per-domain
                approach, resolving cross-domain overlapping hits.

The NAME column is removed from output in all cases; ACC is the canonical ID.
When --iprlookup is set, InterPro and IPR_desc columns are added.

Usage: add_ipr_go.py --tsv FILE --db DBNAME --data-dir DIR [--iprlookup] [--goterms]
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict


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


def load_interpro_names(data_dir):
    """
    Load IPR_ID -> short description from data/interpro_names.tsv.
    Built by install.sh from interpro.xml.gz.
    """
    path = os.path.join(data_dir, "interpro_names.tsv")
    if not os.path.isfile(path):
        return {}
    m = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            tab = line.index("\t")
            m[line[:tab]] = line[tab + 1:]
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


def load_pirsf_thresholds(data_dir):
    """
    Parse data/pirsf/pirsf.dat into PIRSF_ID -> min_bit_score.
    Threshold = max(0, avg_score - 3.5 * std_score), matching InterProScan.
    """
    path = os.path.join(data_dir, "pirsf", "pirsf.dat")
    if not os.path.isfile(path):
        return {}
    m = {}
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            pirsf_id = lines[i][1:].split()[0]  # ID only (before any " child: ..." text)
            if i + 2 < len(lines):
                nums = lines[i + 2].split()
                if len(nums) >= 2:
                    try:
                        avg = float(nums[0])
                        std = float(nums[1])
                        m[pirsf_id] = max(0.0, avg - 3.5 * std)
                    except ValueError:
                        pass
            i += 4
        else:
            i += 1
    return m


def load_pirsf_hierarchy(data_dir):
    """
    Parse pirsf.dat for parent -> set(children) relationships.
    Lines like ">PIRSF000551 child: PIRSF501104 PIRSF501105 ..."
    Returns dict: parent_id -> frozenset(child_ids).
    """
    path = os.path.join(data_dir, "pirsf", "pirsf.dat")
    if not os.path.isfile(path):
        return {}
    hierarchy = {}
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            parts = line.split()
            if len(parts) < 3 or parts[1] != "child:":
                continue
            parent_id = parts[0][1:]
            hierarchy[parent_id] = frozenset(parts[2:])
    return hierarchy


def load_pfam_clans(data_dir):
    """
    Parse pfam_a.dat to build model_name -> clan_id mapping.
    Models that do not belong to a clan are absent from the returned dict.
    """
    path = os.path.join(data_dir, "pfam", "pfam_a.dat")
    if not os.path.isfile(path):
        return {}
    clans = {}
    current_id = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#=GF ID"):
                current_id = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF CL") and current_id:
                clans[current_id] = line.split(None, 2)[2].strip()
            elif line == "//":
                current_id = None
    return clans


def pfam_clan_dedup(rows, clans, overlap_threshold=0.30):
    """
    For each protein, remove Pfam hits that:
      1. Belong to the same clan as an already-accepted (higher-scoring) hit,
         AND
      2. Overlap by more than overlap_threshold of the shorter hit's length.

    Hits without a clan assignment are always kept (never deduplicated against
    each other via clan logic, though they may still be distinct domain
    instances).

    This mirrors InterProScan's Pfam post-processing which resolves redundancy
    within Pfam clans (related families that recognise the same structural fold).
    Rows should be pre-sorted by score descending within each protein so the
    best representative of each clan-region is retained.
    """
    by_protein = defaultdict(list)
    for row in rows:
        by_protein[row.get("target", "")].append(row)

    result = []
    for protein_rows in by_protein.values():
        protein_rows.sort(key=lambda r: float(r.get("score") or 0), reverse=True)
        accepted = []
        for row in protein_rows:
            clan = clans.get(row.get("NAME", ""))
            if clan is None:
                accepted.append(row)
                continue
            try:
                start = int(row.get("target_start") or 0)
                end = int(row.get("target_end") or 0)
            except (ValueError, TypeError):
                accepted.append(row)
                continue
            hit_len = end - start
            redundant = False
            for sel in accepted:
                if clans.get(sel.get("NAME", "")) != clan:
                    continue  # Different clan — no clan-based conflict
                try:
                    s_start = int(sel.get("target_start") or 0)
                    s_end = int(sel.get("target_end") or 0)
                except (ValueError, TypeError):
                    continue
                overlap = max(0, min(end, s_end) - max(start, s_start))
                shorter = min(hit_len, s_end - s_start)
                if shorter > 0 and overlap / shorter > overlap_threshold:
                    redundant = True
                    break
            if not redundant:
                accepted.append(row)
        result.extend(accepted)
    return result


def load_cath_map(data_dir):
    """Load model_id -> G3DSA:x.x.x.x from data/cath/model_to_family_map.tsv."""
    path = os.path.join(data_dir, "cath", "model_to_family_map.tsv")
    if not os.path.isfile(path):
        return {}
    m = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                model_id = parts[0].strip()
                cath_domain = parts[1].strip()
                if model_id and cath_domain:
                    m[model_id] = f"G3DSA:{cath_domain}"
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
    return acc[:dot] if suffix.isdigit() else acc


def cath_base_name(name):
    """Strip -i1/-i2 style suffix: 2vqeL00-i2 -> 2vqeL00"""
    if name and len(name) >= 3 and name[-2] == "i" and name[-3] == "-":
        return name[:-3]
    return name


# ── Filtering helpers ────────────────────────────────────────────────────────

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


def best_hit_per_domain(rows):
    """Keep the highest-scoring hit for each (protein, domain_acc) pair."""
    best = {}
    for row in rows:
        target = row.get("target", "")
        acc = row.get("ACC", "")
        try:
            score = float(row.get("score") or 0)
        except (ValueError, TypeError):
            score = 0.0
        key = (target, acc)
        if key not in best or score > best[key][0]:
            best[key] = (score, row)
    return [row for _, (_, row) in sorted(best.items())]


def non_overlapping_hits(rows, overlap_threshold=0.25):
    """
    From a list of hits (pre-sorted by score descending), greedily select
    those that do not overlap significantly with already-accepted hits.

    Two hits overlap "significantly" when:
      overlap_residues / len(shorter_hit) > overlap_threshold

    This mirrors InterProScan's domain-collapsing behaviour for SUPERFAMILY,
    where multiple HMM models can match the same fold at slightly different
    positions; only non-redundant domain instances are kept.
    """
    accepted = []
    for row in rows:
        try:
            start = int(row.get("target_start") or 0)
            end = int(row.get("target_end") or 0)
        except (ValueError, TypeError):
            accepted.append(row)
            continue
        if end <= start:
            accepted.append(row)
            continue
        hit_len = end - start
        redundant = False
        for sel in accepted:
            try:
                s_start = int(sel.get("target_start") or 0)
                s_end = int(sel.get("target_end") or 0)
            except (ValueError, TypeError):
                continue
            overlap = max(0, min(end, s_end) - max(start, s_start))
            shorter = min(hit_len, s_end - s_start)
            if shorter > 0 and overlap / shorter > overlap_threshold:
                redundant = True
                break
        if not redundant:
            accepted.append(row)
    return accepted


def non_overlapping_per_domain(rows, overlap_threshold=0.25):
    """
    Group rows by (protein, domain_acc), sort each group by score (desc),
    then apply non_overlapping_hits within each group.
    Returns a flat list of kept rows, sorted by (target, score desc).
    """
    groups = defaultdict(list)
    for row in rows:
        key = (row.get("target", ""), row.get("ACC", ""))
        groups[key].append(row)

    result = []
    for group_rows in groups.values():
        group_rows.sort(key=lambda r: float(r.get("score") or 0), reverse=True)
        result.extend(non_overlapping_hits(group_rows, overlap_threshold))
    result.sort(key=lambda r: (r.get("target", ""), -float(r.get("score") or 0)))
    return result


def non_overlapping_per_protein(rows, overlap_threshold=0.35):
    """
    Apply non-overlapping greedy selection globally across ALL hits for each
    protein, regardless of which domain/family they belong to.

    This matches InterProScan's SUPERFAMILY post-processing (ass3.pl,
    $percentsame = 0.35): hits from *different* superfamilies that map to
    the same protein region are collapsed — only the highest-scoring hit per
    region survives, irrespective of which SSF family it represents.
    """
    by_protein = defaultdict(list)
    for row in rows:
        by_protein[row.get("target", "")].append(row)

    result = []
    for protein_rows in by_protein.values():
        protein_rows.sort(key=lambda r: float(r.get("score") or 0), reverse=True)
        result.extend(non_overlapping_hits(protein_rows, overlap_threshold))
    result.sort(key=lambda r: (r.get("target", ""), -float(r.get("score") or 0)))
    return result


def load_sfld_hierarchy(data_dir):
    """
    Parse data/sfld/sfld_hierarchy.tsv into a mapping:
      child_id -> set(ancestor_ids)

    The TSV has up to 3 columns representing a path from root to leaf:
      SFLDS (superfamily)  SFLDG (group)  [SFLDF (family)]

    Each unique ID that appears in a non-leftmost column is considered a
    child of everything to its left on the same row.  We collect all
    ancestor IDs for each child so the filter can simply ask
    "does this ID have any descendant also in the matched set?"

    Returns: (children_map, ancestors_map)
      children_map:  parent_id -> frozenset(direct_child_ids)
      ancestors_map: child_id  -> frozenset(all_ancestor_ids)
    """
    path = os.path.join(data_dir, "sfld", "sfld_hierarchy.tsv")
    if not os.path.isfile(path):
        return {}, {}

    children_map = defaultdict(set)
    ancestors_map = defaultdict(set)

    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            cols = [c.strip() for c in cols if c.strip()]
            for i in range(1, len(cols)):
                parent = cols[i - 1]
                child = cols[i]
                children_map[parent].add(child)
                # child inherits all ancestors of parent plus parent itself
                ancestors_map[child].add(parent)
                ancestors_map[child].update(ancestors_map.get(parent, set()))

    return (
        {k: frozenset(v) for k, v in children_map.items()},
        {k: frozenset(v) for k, v in ancestors_map.items()},
    )


def filter_sfld_hierarchy(rows, ancestors_map):
    """
    For each protein, remove SFLD entries that are ancestors of a more
    specific entry that also matched.

    Example: if a protein matches SFLDS00014 (superfamily), SFLDG00301
    (group), and SFLDF00157 (family), only SFLDF00157 is kept because it is
    the most specific and the others are its ancestors.
    """
    if not ancestors_map:
        return rows
    by_protein = defaultdict(list)
    for row in rows:
        by_protein[row.get("target", "")].append(row)
    result = []
    for protein_rows in by_protein.values():
        matched_accs = {r.get("ACC", "") for r in protein_rows}
        kept = []
        for row in protein_rows:
            acc = row.get("ACC", "")
            # Keep only if none of its descendants are also matched
            has_matched_descendant = any(
                acc in ancestors_map.get(other, frozenset())
                for other in matched_accs if other != acc
            )
            if not has_matched_descendant:
                kept.append(row)
        result.extend(kept)
    return result


def filter_pirsf_hierarchy(rows, hierarchy):
    """
    For each protein, remove PIRSF entries whose direct children also matched
    (keep only the most specific family in the PIRSF tree).
    """
    if not hierarchy:
        return rows
    by_protein = defaultdict(list)
    for row in rows:
        by_protein[row.get("target", "")].append(row)
    result = []
    for protein_rows in by_protein.values():
        matched = {r.get("ACC", "") for r in protein_rows}
        kept = [r for r in protein_rows
                if not (hierarchy.get(r.get("ACC", ""), frozenset()) & matched)]
        result.extend(kept)
    return result


# ── IPR / GO / IPR_desc annotation ───────────────────────────────────────────

def annotate_ipr_go(rows, db_xml, sig2ipr, ipr2go, ipr_names, do_ipr, do_go):
    for row in rows:
        acc = (row.get("ACC") or "").strip()
        ipr_list = sig2ipr.get((db_xml, acc), [])
        if not ipr_list and acc and "." in acc:
            ipr_list = sig2ipr.get((db_xml, acc.split(".")[0]), [])
        if do_ipr:
            row["InterPro"] = ",".join(ipr_list) if ipr_list else ""
            if ipr_names is not None:
                row["IPR_desc"] = ",".join(
                    ipr_names[ipr] for ipr in ipr_list if ipr in ipr_names
                )
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

# Databases with custom post-processing (always routed through this script)
CUSTOM_DBS = {"panther", "hamap", "superfamily", "pfam", "ncbifam",
              "sfld", "pirsf", "cath"}


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

    needs_custom = db in CUSTOM_DBS

    if not args.iprlookup and not args.goterms and not needs_custom:
        with open(args.tsv) as f:
            sys.stdout.write(f.read())
        return

    # Load IPR/GO/names data if needed
    sig2ipr, ipr2go, ipr_names = {}, {}, None
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
        if args.iprlookup:
            ipr_names = load_interpro_names(data_dir)

    with open(args.tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames)
        if args.iprlookup and "InterPro" not in fieldnames:
            fieldnames.append("InterPro")
            fieldnames.append("IPR_desc")
        elif args.iprlookup and "IPR_desc" not in fieldnames:
            idx = fieldnames.index("InterPro")
            fieldnames.insert(idx + 1, "IPR_desc")
        if args.goterms and "GO" not in fieldnames:
            fieldnames.append("GO")
        # Remove NAME column from output - ACC is the canonical identifier
        fieldnames = [col for col in fieldnames if col != "NAME"]

        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")

        # ── PANTHER ──────────────────────────────────────────────────────────
        if db == "panther":
            panther_names = load_panther_names(data_dir)
            all_rows = list(reader)
            for row in all_rows:
                base_id = pthr_base_id(row.get("NAME") or row.get("ACC") or "")
                row["ACC"] = base_id
                if "DESC" in fieldnames and not (row.get("DESC") or "").strip():
                    row["DESC"] = panther_names.get(base_id, "")
            rows = best_hit_per_target(all_rows)
            annotate_ipr_go(rows, db_xml, sig2ipr, ipr2go, ipr_names, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(rows)
            return

        # ── HAMAP ─────────────────────────────────────────────────────────────
        if db == "hamap":
            all_rows = list(reader)
            for row in all_rows:
                if not (row.get("ACC") or "").strip():
                    row["ACC"] = (row.get("NAME") or "").strip()
            rows = best_hit_per_target(all_rows)
            annotate_ipr_go(rows, db_xml, sig2ipr, ipr2go, ipr_names, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(rows)
            return

        # ── SUPERFAMILY ───────────────────────────────────────────────────────
        # InterProScan (ass3.pl) uses E < 1e-4 as the stated default cutoff,
        # but the *effective* cutoff is adjusted per-protein by:
        #   tempcutoff = cutoff / (n_sfs / totmods)
        # where n_sfs is the number of unique superfamilies hit by that protein
        # and totmods is the total number of HMM models in the database.  For a
        # typical protein matching only a few dozen superfamilies out of ~2,000+
        # SUPERFAMILY HMM models, this ratio is <<1, making the effective cutoff
        # 10–100× more permissive than 1e-4.  Using 1e-3 here is a reasonable
        # approximation of that adjusted threshold and helps recall.
        #
        # Overlap filter: global per-protein non-overlapping selection with
        # overlap_threshold = 0.35 ($percentsame = 0.35 in ass3.pl).
        # Applied across ALL SSF hits for a protein (not per-domain), so
        # cross-family overlapping hits are also collapsed.
        if db == "superfamily":
            all_rows = list(reader)
            filtered = []
            for row in all_rows:
                try:
                    evalue = float(row.get("evalue") or 1)
                except (ValueError, TypeError):
                    evalue = 1.0
                if evalue > 1e-3:
                    continue
                clean = ssf_id(row.get("ACC") or row.get("NAME") or "")
                row["ACC"] = clean
                filtered.append(row)
            rows = non_overlapping_per_protein(filtered, overlap_threshold=0.35)
            annotate_ipr_go(rows, db_xml, sig2ipr, ipr2go, ipr_names, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(rows)
            return

        # ── Pfam ──────────────────────────────────────────────────────────────
        # Step 1: domain-level GA threshold (per-hit score filter).
        # Step 2: sequence-level GA threshold — for each (protein, model) pair,
        #   sum all domain scores as a proxy for HMMER's sequence-level score;
        #   drop all domains for that pair if the sum is below GA_seq.
        #   (GA_seq and GA_dom are equal for most Pfam models; this catches the
        #    minority where a single weak domain passes GA_dom but the protein
        #    lacks enough cumulative evidence.)
        # Step 3: clan-based overlap deduplication — within each protein, if two
        #   hits from the same Pfam clan overlap by > 30% of the shorter hit,
        #   only the higher-scoring one is kept.  This matches InterProScan's
        #   Pfam post-processing which resolves redundancy within clan families.
        if db == "pfam":
            ga = load_thresholds(os.path.join(data_dir, "pfam", "pfam_a.ga"))
            clans = load_pfam_clans(data_dir)
            all_rows = list(reader)

            # Step 1: domain-level GA
            after_dom = []
            for row in all_rows:
                name = row.get("NAME", "")
                try:
                    score = float(row.get("score") or 0)
                except (ValueError, TypeError):
                    score = 0.0
                if name in ga and score < ga[name][1]:
                    continue
                row["ACC"] = strip_version(row.get("ACC") or "")
                after_dom.append(row)

            # Step 2: sequence-level GA (sum domain scores per protein+model)
            if ga:
                by_protmod = defaultdict(list)
                for row in after_dom:
                    by_protmod[(row.get("target", ""), row.get("NAME", ""))].append(row)
                after_seq = []
                for (_, name), group in by_protmod.items():
                    seq_thresh = ga.get(name, (None, None))[0]
                    if seq_thresh is not None:
                        total = sum(float(r.get("score") or 0) for r in group)
                        if total < seq_thresh:
                            continue
                    after_seq.extend(group)
            else:
                after_seq = after_dom

            # Step 3: clan-based overlap dedup.
            # Threshold 0.50: only remove a clan member if it overlaps the
            # best-scoring member of the same clan by more than half its length.
            # A looser threshold (vs. the earlier 0.30) avoids discarding
            # adjacent same-clan domains that genuinely co-occur in tandem
            # repeats or multi-domain proteins, which was causing false negatives.
            kept = pfam_clan_dedup(after_seq, clans, overlap_threshold=0.50)

            annotate_ipr_go(kept, db_xml, sig2ipr, ipr2go, ipr_names, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(kept)
            return

        # ── NCBIfam ───────────────────────────────────────────────────────────
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
                if name in tc and score < tc[name][0]:
                    continue
                raw_acc = row.get("ACC") or row.get("NAME") or ""
                row["ACC"] = strip_version(raw_acc)
                kept.append(row)
            annotate_ipr_go(kept, db_xml, sig2ipr, ipr2go, ipr_names, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(kept)
            return

        # ── SFLD ──────────────────────────────────────────────────────────────
        # InterProScan's sfld_postprocess binary additionally performs
        # residue-level active-site matching (using sfld_sites.annot and
        # HMMER's Stockholm alignment output) to discard domain hits that
        # pass GA thresholds but lack the catalytic residues.  Nail does not
        # produce alignment files in a compatible format, so full site-level
        # filtering is not available here.
        #
        # We apply two improvements that are feasible without alignment data:
        #   1. GA domain-level threshold (same as before).
        #   2. Hierarchy filter: SFLD is organised as SFLDS > SFLDG > SFLDF.
        #      When a protein matches both a parent and a child entry, only the
        #      most specific (child) match is retained — matching the same
        #      logic used for PIRSF.  Data from sfld_hierarchy.tsv shipped
        #      with InterProScan.
        if db == "sfld":
            ga = load_thresholds(os.path.join(data_dir, "sfld", "sfld.ga"))
            _, ancestors_map = load_sfld_hierarchy(data_dir)
            all_rows = list(reader)
            kept = []
            for row in all_rows:
                name = row.get("NAME", "")
                try:
                    score = float(row.get("score") or 0)
                except (ValueError, TypeError):
                    score = 0.0
                if name in ga and score < ga[name][1]:
                    continue
                kept.append(row)
            kept = filter_sfld_hierarchy(kept, ancestors_map)
            annotate_ipr_go(kept, db_xml, sig2ipr, ipr2go, ipr_names, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(kept)
            return

        # ── PIRSF ─────────────────────────────────────────────────────────────
        # Two-stage filtering:
        #   1. Per-family score threshold from pirsf.dat (avg - 3.5 * std).
        #   2. Hierarchy filter: when a protein matches both a parent and a child
        #      PIRSF family, only the child (more specific) is reported.
        # Note: InterProScan's full hierarchy-aware algorithm also considers
        # protein length and HMM coverage; this is a best-effort approximation.
        if db == "pirsf":
            thresh_file = os.path.join(data_dir, "pirsf", "pirsf.thresh")
            thresh_tuples = load_thresholds(thresh_file)
            if thresh_tuples:
                thresholds = {k: v[0] for k, v in thresh_tuples.items()}
            else:
                thresholds = load_pirsf_thresholds(data_dir)
            hierarchy = load_pirsf_hierarchy(data_dir)
            all_rows = list(reader)
            kept = []
            for row in all_rows:
                acc = (row.get("ACC") or "").strip()
                try:
                    score = float(row.get("score") or 0)
                except (ValueError, TypeError):
                    score = 0.0
                if acc in thresholds and score < thresholds[acc]:
                    continue
                if not (row.get("DESC") or "").strip():
                    row["DESC"] = (row.get("NAME") or "").strip()
                kept.append(row)
            kept = filter_pirsf_hierarchy(kept, hierarchy)
            annotate_ipr_go(kept, db_xml, sig2ipr, ipr2go, ipr_names, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(kept)
            return

        # ── CATH (Gene3D) ─────────────────────────────────────────────────────
        # InterProScan uses the `cath-resolve-hits` binary which finds the
        # non-overlapping set of domain hits that *maximises the total bit
        # score* via dynamic programming (not greedy best-per-domain).
        # Default filters applied by cath-resolve-hits:
        #   --worst-permissible-evalue 0.001
        #   --worst-permissible-bitscore 10
        #   --overlap-trim-spec 30/10  (allows ~10-residue boundary overlap)
        #
        # Our previous approach (best_hit_per_domain) only collapsed redundancy
        # within the same G3DSA domain but left cross-domain overlapping hits
        # intact, causing false positives when two different G3DSA families
        # matched the same protein region.
        #
        # Replacement: global per-protein non-overlapping selection
        # (same algorithm as SUPERFAMILY post-processing), which resolves
        # cross-domain overlaps by keeping the higher-scoring hit per region.
        # Overlap threshold 0.20 is stricter than cath-resolve-hits' effective
        # ~10% trim, but appropriate given that our greedy approach cannot
        # optimise globally like the DP solver.
        if db == "cath":
            cath_map = load_cath_map(data_dir)
            all_rows = list(reader)
            mapped = []
            for row in all_rows:
                raw_name = row.get("ACC") or row.get("NAME") or ""
                base = cath_base_name(raw_name)
                domain = cath_map.get(base, "")
                if not domain:
                    continue
                try:
                    evalue = float(row.get("evalue") or 1)
                    score = float(row.get("score") or 0)
                except (ValueError, TypeError):
                    evalue = 1.0
                    score = 0.0
                if evalue > 1e-4:
                    continue
                if score < 10:  # cath-resolve-hits --worst-permissible-bitscore 10
                    continue
                row["ACC"] = domain
                mapped.append(row)
            rows = non_overlapping_per_protein(mapped, overlap_threshold=0.20)
            annotate_ipr_go(rows, db_xml, sig2ipr, ipr2go, ipr_names, args.iprlookup, args.goterms)
            writer.writeheader()
            writer.writerows(rows)
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
                row["IPR_desc"] = ",".join(
                    ipr_names[ipr] for ipr in ipr_list
                    if ipr_names and ipr in ipr_names
                ) if ipr_names else ""
            if args.goterms:
                go_set = set()
                for ipr in ipr_list:
                    go_set.update(ipr2go.get(ipr, []))
                row["GO"] = ",".join(sorted(go_set)) if go_set else ""
            writer.writerow(row)


if __name__ == "__main__":
    main()
