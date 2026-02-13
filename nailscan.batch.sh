#!/usr/bin/env bash
# NailScan batch: combine all .fa/.faa/.fasta in a directory, run nail Pfam search, split to results/pfam/{genome_id}.tsv.gz
# Usage: ./nailscan.batch.sh [-f fasta_dir] [-h pfam_hmm] [-o results_dir] [-t N]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAIL_THREADS="${NAIL_THREADS:-$(nproc 2>/dev/null || echo 16)}"
PFAM_HMM="${SCRIPT_DIR}/data/pfam/pfam_a.hmm"
FASTA_DIR="${SCRIPT_DIR}/fasta"
RESULTS_PFAM="${SCRIPT_DIR}/results/pfam"

usage() {
  echo "Usage: $0 [-f fasta_dir] [-h pfam_hmm] [-o results_dir] [-t N]"
  echo "  -f  Directory containing .fa/.faa/.fasta files (default: fasta/)"
  echo "  -h  Pfam HMM path (default: data/pfam/pfam_a.hmm)"
  echo "  -o  Output directory for per-genome .tsv.gz (default: results/pfam/)"
  echo "  -t  Threads (default: nproc or 16)"
}
[[ " $* " == *" --help "* ]] || [[ " $* " == *" -help "* ]] && { usage; exit 0; }

# Parse -h (HMM), -f (fasta dir), -o (output dir), -t (threads)
while getopts "h:f:o:t:" opt; do
  case $opt in
    h) PFAM_HMM="$OPTARG" ;;
    f) FASTA_DIR="$OPTARG" ;;
    o) RESULTS_PFAM="$OPTARG" ;;
    t) NAIL_THREADS="$OPTARG" ;;
    *) usage; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

# Validate Pfam HMM file exists
if [[ ! -f "$PFAM_HMM" ]]; then
  echo "Error: Pfam HMM not found: $PFAM_HMM"
  exit 1
fi
# Validate input directory exists
if [[ ! -d "$FASTA_DIR" ]]; then
  echo "Error: Fasta directory not found: $FASTA_DIR"
  exit 1
fi
# Ensure nail is on PATH (conda env)
if ! command -v nail &>/dev/null; then
  echo "Error: nail not in PATH. Activate env: conda activate nail"
  exit 1
fi

mkdir -p "$RESULTS_PFAM"
PFAM_HMM="$(cd "$(dirname "$PFAM_HMM")" && pwd)/$(basename "$PFAM_HMM")"
FASTA_DIR="$(cd "$FASTA_DIR" && pwd)"
RESULTS_PFAM="$(cd "$RESULTS_PFAM" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

BATCH_FASTA="${WORK}/batch.fa"
MAP_FILE="${PFAM_HMM}.map"

# Collect all .fa, .faa, .fasta files (each file once) and build combined fasta with headers genome_id|seqid
printf 'building combined fasta... '
t_start=$(date +%s.%N)
: > "$BATCH_FASTA"
seen=""
for fa in "$FASTA_DIR"/*.fa "$FASTA_DIR"/*.faa "$FASTA_DIR"/*.fasta; do
  [[ -f "$fa" ]] || continue
  [[ " $seen " == *" $fa "* ]] && continue
  seen="$seen $fa"
  base="$(basename "$fa")"
  genome_id="$base"
  [[ "$base" == *.fasta ]] && genome_id="${base%.fasta}"
  [[ "$base" == *.faa ]] && genome_id="${base%.faa}"
  [[ "$base" == *.fa ]] && genome_id="${base%.fa}"
  awk -v g="$genome_id" '/^>/ { print ">" g "|" substr($0,2); next } { print }' "$fa" >> "$BATCH_FASTA"
done
n_seqs=$(grep -c '^>' "$BATCH_FASTA" || true)
t_end=$(date +%s.%N)
printf "done (%.2fs)\n" "$(awk -v s="$t_start" -v e="$t_end" 'BEGIN { printf "%.2f", e-s }')"
echo "    ${n_seqs} sequences"

# Run nail search (HMM vs combined fasta) with thread count
(cd "$WORK" && nail search -t "$NAIL_THREADS" --tbl-out results.tbl "$PFAM_HMM" "$BATCH_FASTA")

if [[ ! -f "$WORK/results.tbl" ]]; then
  echo "Error: nail did not produce results.tbl"
  exit 1
fi

# Build Pfam NAME -> ACC,DESC map from HMM if not present
if [[ ! -f "$MAP_FILE" ]]; then
  printf 'building pfam map...        '
  t_start=$(date +%s.%N)
  if [[ -f "${SCRIPT_DIR}/scripts/make_pfam_map.sh" ]]; then
    bash "${SCRIPT_DIR}/scripts/make_pfam_map.sh" "$PFAM_HMM" "$MAP_FILE" >/dev/null
  else
    awk 'BEGIN{name="";acc="";desc=""}
         /^NAME[ \t]+/ { if (name != "") print name "\t" acc "\t" desc; name=$2; acc=""; desc="" }
         /^ACC[ \t]+/  { acc=$2 }
         /^DESC[ \t]+/ { desc=$0; sub(/^DESC[ \t]+/, "", desc) }
         END { if (name != "") print name "\t" acc "\t" desc }' "$PFAM_HMM" > "$MAP_FILE"
  fi
  t_end=$(date +%s.%N)
  printf "done (%.2fs)\n" "$(awk -v s="$t_start" -v e="$t_end" 'BEGIN { printf "%.2f", e-s }')"
fi

HEADER_NO_MAP="target	NAME	target_start	target_end	query_start	query_end	score	bias	evalue	cell_frac"
HEADER_MAP="target	NAME	ACC	DESC	target_start	target_end	query_start	query_end	score	bias	evalue	cell_frac"
USE_MAP=0
[[ -f "$MAP_FILE" ]] && USE_MAP=1
mkdir -p "${WORK}/per_genome"

# Split combined results by genome_id (prefix before |) into per-genome TSVs, with optional ACC/DESC
printf 'splitting results...        '
t_start=$(date +%s.%N)
if [[ $USE_MAP -eq 1 ]]; then
  awk -v work="${WORK}/per_genome" -v hdr="$HEADER_MAP" '
    NR==FNR { map[$1]=$2"\t"$3; next }
    /^#/ || NF<10 { next }
    {
      target=$1; query=$2
      idx=index(target,"|"); if (idx==0) next
      genome_id=substr(target,1,idx-1)
      gene_id=substr(target,idx+1)
      f=work "/" genome_id ".tsv"
      if (!(genome_id in seen)) { print hdr > f; seen[genome_id]=1 }
      acc_desc=map[query]; if (acc_desc=="") acc_desc="\t"
      split(acc_desc,a,"\t"); acc=a[1]; desc=a[2]
      print gene_id"\t"query"\t"acc"\t"desc"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 >> f
    }
  ' "$MAP_FILE" "$WORK/results.tbl"
else
  awk -v work="${WORK}/per_genome" -v hdr="$HEADER_NO_MAP" '
    /^#/ || NF<10 { next }
    {
      target=$1
      idx=index(target,"|"); if (idx==0) next
      genome_id=substr(target,1,idx-1)
      gene_id=substr(target,idx+1)
      f=work "/" genome_id ".tsv"
      if (!(genome_id in seen)) { print hdr > f; seen[genome_id]=1 }
      print gene_id"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 >> f
    }
  ' "$WORK/results.tbl"
fi
# Gzip each per-genome TSV into results directory
for tsv in "${WORK}/per_genome"/*.tsv; do
  [[ -f "$tsv" ]] || continue
  gzip -c "$tsv" > "${RESULTS_PFAM}/$(basename "$tsv" .tsv).tsv.gz"
done
t_end=$(date +%s.%N)
printf "done (%.2fs)\n" "$(awk -v s="$t_start" -v e="$t_end" 'BEGIN { printf "%.2f", e-s }')"

echo ""
echo "Results in ${RESULTS_PFAM}/"
