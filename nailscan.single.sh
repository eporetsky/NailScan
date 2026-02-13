#!/usr/bin/env bash
# NailScan single: run nail Pfam search on one FASTA file; write TSV with ACC/DESC when map exists
# Usage: ./nailscan.single.sh -f <input.fasta> [-h pfam_hmm] [-o output_dir] [-t N]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAIL_THREADS="${NAIL_THREADS:-$(nproc 2>/dev/null || echo 16)}"
PFAM_HMM="${SCRIPT_DIR}/data/pfam/pfam_a.hmm"
INPUT_FASTA=""
OUTPUT_DIR="${SCRIPT_DIR}/results"

usage() {
  echo "Usage: $0 -f <input.fasta> [-h pfam_hmm] [-o output_dir] [-t N]"
  echo "  -f  FASTA file (required)"
  echo "  -h  Pfam HMM path (default: data/pfam/pfam_a.hmm)"
  echo "  -o  Output directory (default: results/)"
  echo "  -t  Threads (default: nproc or 16)"
}
[[ " $* " == *" --help "* ]] || [[ " $* " == *" -help "* ]] && { usage; exit 0; }

# Parse -h (HMM), -f (fasta), -o (output), -t (threads)
while getopts "h:f:o:t:" opt; do
  case $opt in
    h) PFAM_HMM="$OPTARG" ;;
    f) INPUT_FASTA="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    t) NAIL_THREADS="$OPTARG" ;;
    *) usage; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

if [[ -z "$INPUT_FASTA" ]]; then
  echo "Error: -f <input.fasta> is required"
  usage
  exit 1
fi
if [[ ! -f "$PFAM_HMM" ]]; then
  echo "Error: Pfam HMM not found: $PFAM_HMM"
  exit 1
fi
if [[ ! -f "$INPUT_FASTA" ]]; then
  echo "Error: Fasta file not found: $INPUT_FASTA"
  exit 1
fi
if ! command -v nail &>/dev/null; then
  echo "Error: nail not in PATH. Activate env: conda activate nail"
  exit 1
fi

# Resolve absolute paths; output file named from input basename
PFAM_HMM="$(cd "$(dirname "$PFAM_HMM")" && pwd)/$(basename "$PFAM_HMM")"
INPUT_FASTA="$(cd "$(dirname "$INPUT_FASTA")" && pwd)/$(basename "$INPUT_FASTA")"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
base="$(basename "$INPUT_FASTA")"
genome_id="${base%.*}"
[[ "$base" == *.fasta ]] && genome_id="${base%.fasta}"
[[ "$base" == *.faa ]] && genome_id="${base%.faa}"
[[ "$base" == *.fa ]] && genome_id="${base%.fa}"
OUT_TSV="${OUTPUT_DIR}/${genome_id}.pfam.tsv"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Run nail search (HMM vs single fasta) with requested thread count
echo "Running: $INPUT_FASTA (threads: $NAIL_THREADS)"
(cd "$WORK" && nail search -t "$NAIL_THREADS" --tbl-out results.tbl "$PFAM_HMM" "$INPUT_FASTA")

if [[ ! -f "$WORK/results.tbl" ]]; then
  echo "Error: nail did not produce results.tbl"
  exit 1
fi

MAP_FILE="${PFAM_HMM}.map"
# Build Pfam NAME -> ACC,DESC map from HMM file if not already present
if [[ ! -f "$MAP_FILE" ]]; then
  if [[ -f "${SCRIPT_DIR}/scripts/make_pfam_map.sh" ]]; then
    bash "${SCRIPT_DIR}/scripts/make_pfam_map.sh" "$PFAM_HMM" "$MAP_FILE" >/dev/null
  else
    awk 'BEGIN{name="";acc="";desc=""}
         /^NAME[ \t]+/ { if (name != "") print name "\t" acc "\t" desc; name=$2; acc=""; desc="" }
         /^ACC[ \t]+/  { acc=$2 }
         /^DESC[ \t]+/ { desc=$0; sub(/^DESC[ \t]+/, "", desc) }
         END { if (name != "") print name "\t" acc "\t" desc }' "$PFAM_HMM" > "$MAP_FILE"
  fi
fi

HEADER_NO_MAP="target	NAME	target_start	target_end	query_start	query_end	score	bias	evalue	cell_frac"
HEADER_MAP="target	NAME	ACC	DESC	target_start	target_end	query_start	query_end	score	bias	evalue	cell_frac"
# Write TSV with header; add ACC and DESC from map when available
if [[ -f "$MAP_FILE" ]]; then
  awk -v hdr="$HEADER_MAP" '
    NR==FNR { map[$1]=$2"\t"$3; next }
    BEGIN { print hdr }
    /^#/ || NF<10 { next }
    {
      query=$2
      acc_desc=map[query]; if (acc_desc=="") acc_desc="\t"
      split(acc_desc,a,"\t"); acc=a[1]; desc=a[2]
      print $1"\t"query"\t"acc"\t"desc"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10
    }
  ' "$MAP_FILE" "$WORK/results.tbl" > "$OUT_TSV"
else
  awk -v hdr="$HEADER_NO_MAP" '
    BEGIN { print hdr }
    /^#/ || NF<10 { next }
    { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }
  ' "$WORK/results.tbl" > "$OUT_TSV"
fi

echo "Wrote: $OUT_TSV"
echo "Done. Results: $OUT_TSV"
