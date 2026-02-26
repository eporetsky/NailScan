#!/usr/bin/env bash
# NailScan batch wrapper — annotate multiple genome FASTAs by delegating each
# batch to nailscan.single.sh, then splitting the combined output into
# per-genome TSV.gz files.
#
# All HMM searching, filtering, and column manipulation is handled by
# nailscan.single.sh; this script only handles:
#   1. Chunking a directory of FASTAs into batches (amortises the cost of
#      reading large HMM databases — one read per batch, not per genome).
#   2. Prefixing protein IDs with the genome name so a single nail run covers
#      an entire batch.
#   3. Splitting the combined output back into per-genome TSV.gz files.
#
# Usage:
#   ./nailscan.batch.sh -b <base> [-f fasta_dir] [-o results_dir]
#                       [--batch-size N] [-t threads] [-d data_dir]
#                       [--appl db1,db2] [--iprlookup] [--goterms]
#
# --batch-size controls the RAM / query-DB-read tradeoff:
#   0 (default): all genomes in one call — HMM database read once, highest RAM
#   100:         100 genomes per call — HMM database read ceil(N/100) times
#   1:           one genome per call — equivalent to looping over nailscan.single.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SINGLE_SH="${SCRIPT_DIR}/nailscan.single.sh"

# ── Defaults ──────────────────────────────────────────────────────────────────
FASTA_DIR="${SCRIPT_DIR}/fasta"
RESULTS_ROOT="${SCRIPT_DIR}/results"
OUTPUT_FILE_BASE=""
BATCH_SIZE=0        # 0 = all genomes in one call
PASSTHROUGH=()      # flags forwarded unchanged to nailscan.single.sh

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  echo "Usage: $0 -b <output_base> [-f fasta_dir] [-o results_dir] [--batch-size N]"
  echo "          [-t threads] [-d data_dir] [--appl db1,db2] [--iprlookup] [--goterms]"
  echo ""
  echo "  -b / --output-file-base   Output subdir name under <results_dir>/ (required)"
  echo "  -f                        Directory with .fa/.faa/.fasta files (default: fasta/)"
  echo "  -o                        Root results directory (default: results/)"
  echo "  --batch-size N            Genomes per nailscan.single.sh call (default: 0 = all)."
  echo "                            Lower values reduce peak RAM at the cost of re-reading"
  echo "                            large HMM databases more often."
  echo "  -t / -d / --appl / --iprlookup / --goterms  Passed through to nailscan.single.sh"
  echo "  --help                    This message"
}

[[ " $* " == *" --help "* || " $* " == *" -help "* ]] && { usage; exit 0; }

# ── Arg parsing ───────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--output-file-base) OUTPUT_FILE_BASE="$2"; shift 2 ;;
    -f)                    FASTA_DIR="$2"; shift 2 ;;
    -o)                    RESULTS_ROOT="$2"; shift 2 ;;
    --batch-size)          BATCH_SIZE="$2"; shift 2 ;;
    # Forward these directly to nailscan.single.sh
    -t|--threads)          PASSTHROUGH+=("-t" "$2"); shift 2 ;;
    -d|--data-dir)         PASSTHROUGH+=("-d" "$2"); shift 2 ;;
    --appl|--applications) PASSTHROUGH+=("--appl" "$2"); shift 2 ;;
    --iprlookup|--goterms) PASSTHROUGH+=("$1"); shift ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

# ── Validation ────────────────────────────────────────────────────────────────
if [[ -z "$OUTPUT_FILE_BASE" ]]; then
  echo "Error: -b <output_file_base> is required"; usage; exit 1
fi
if [[ ! -f "$SINGLE_SH" ]]; then
  echo "Error: nailscan.single.sh not found at $SINGLE_SH"; exit 1
fi
if [[ ! -d "$FASTA_DIR" ]]; then
  echo "Error: FASTA directory not found: $FASTA_DIR"; exit 1
fi

FASTA_DIR="$(cd "$FASTA_DIR" && pwd)"
mkdir -p "$RESULTS_ROOT/$OUTPUT_FILE_BASE"
OUTPUT_DIR="$(cd "$RESULTS_ROOT/$OUTPUT_FILE_BASE" && pwd)"

# ── Collect FASTA files ───────────────────────────────────────────────────────
FA_LIST=()
for fa in "$FASTA_DIR"/*.fa "$FASTA_DIR"/*.faa "$FASTA_DIR"/*.fasta; do
  [[ -f "$fa" ]] && FA_LIST+=("$fa")
done
n_genomes=${#FA_LIST[@]}
if [[ $n_genomes -eq 0 ]]; then
  echo "No FASTA files (.fa / .faa / .fasta) found in ${FASTA_DIR}"; exit 1
fi

# ── Determine chunk size ──────────────────────────────────────────────────────
chunk="${BATCH_SIZE:-0}"
if [[ "$chunk" -le 0 ]] || [[ "$chunk" -ge "$n_genomes" ]]; then
  chunk="$n_genomes"
fi
n_chunks=$(( (n_genomes + chunk - 1) / chunk ))
if [[ "$n_chunks" -gt 1 ]]; then
  echo "Processing ${n_genomes} genome(s) in ${n_chunks} batch(es) of up to ${chunk} (--batch-size ${chunk})"
else
  echo "Processing ${n_genomes} genome(s) in a single batch"
fi

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# ── Main loop ─────────────────────────────────────────────────────────────────
batch_n=0
batch_start=0
while [[ $batch_start -lt $n_genomes ]]; do
  batch_n=$(( batch_n + 1 ))
  batch_end=$(( batch_start + chunk ))
  [[ $batch_end -gt $n_genomes ]] && batch_end=$n_genomes

  echo ""
  echo "=== Batch ${batch_n}/${n_chunks} (genomes $((batch_start+1))-${batch_end}) ==="

  # Build a combined FASTA: prefix every header with "genome_id|"
  BATCH_FASTA="${WORK}/batch_${batch_n}.fa"
  : > "$BATCH_FASTA"
  for (( i=batch_start; i<batch_end; i++ )); do
    fa="${FA_LIST[$i]}"
    genome_id="$(basename "$fa" | sed 's/\.[^.]*$//')"
    awk -v g="$genome_id" '/^>/{print ">" g "|" substr($0,2); next} {print}' "$fa"
  done >> "$BATCH_FASTA"

  # Run nailscan.single.sh on the combined FASTA
  BATCH_OUTDIR="${WORK}/out_${batch_n}"
  mkdir -p "$BATCH_OUTDIR"

  "$SINGLE_SH" \
    -f "$BATCH_FASTA" \
    -b "batch_${batch_n}" \
    -o "$BATCH_OUTDIR" \
    "${PASSTHROUGH[@]+"${PASSTHROUGH[@]}"}"

  BATCH_TSV="${BATCH_OUTDIR}/batch_${batch_n}.tsv"
  if [[ ! -f "$BATCH_TSV" ]]; then
    echo "Warning: no output produced for batch ${batch_n}, skipping"
    batch_start=$batch_end
    continue
  fi

  # Split the combined TSV back into per-genome TSVs.
  # The target column contains "genome_id|gene_id"; strip the prefix and route
  # each row to the matching output file. Rows for new genomes get a header.
  python3 - "$BATCH_TSV" "$OUTPUT_DIR" <<'PYEOF'
import csv, os, sys
src, out_dir = sys.argv[1], sys.argv[2]
files = {}  # genome_id -> (file_handle, csv_writer)
with open(src, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    headers = list(reader.fieldnames or [])
    for row in reader:
        target = row.get("target", "")
        if "|" not in target:
            continue
        genome_id, gene_id = target.split("|", 1)
        row["target"] = gene_id
        if genome_id not in files:
            path = os.path.join(out_dir, f"{genome_id}.tsv")
            first_write = not os.path.exists(path) or os.path.getsize(path) == 0
            fh = open(path, "a", newline="")
            w = csv.DictWriter(fh, fieldnames=headers, delimiter="\t",
                               extrasaction="ignore", lineterminator="\n")
            if first_write:
                w.writeheader()
            files[genome_id] = (fh, w)
        files[genome_id][1].writerow(row)
for fh, _ in files.values():
    fh.close()
PYEOF

  batch_start=$batch_end
done

# Compress all per-genome TSVs produced in this run
shopt -s nullglob
for tsv in "$OUTPUT_DIR"/*.tsv; do
  gzip -f "$tsv"
done
shopt -u nullglob

echo ""
echo "Done. Results: ${OUTPUT_DIR}/"
