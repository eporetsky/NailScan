#!/usr/bin/env bash
# NailScan batch: combine FASTAs, run nail for each selected DB, split to per-genome TSV.gz.
# By default run all DBs from config. Use --applications to limit.
# Usage: ./nailscan.batch.sh -b <output_file_base> [-f fasta_dir] [-o results_dir] [-t N] [--batch-size N] [--applications ...] [--iprlookup] [--goterms]
# Per-DB outputs are temporary; merged TSV per genome: <results_dir>/<output_file_base>/<genome_id>.tsv.gz
#
# --batch-size controls the RAM / query-DB-read tradeoff:
#   default (0 = unlimited): all genomes in one nail call — query DB read once, highest RAM
#   --batch-size 100:        100 genomes per nail call — query DB read ceil(N/100) times, lower RAM
#   --batch-size 1:          equivalent to nailscan.single.sh in a loop — query DB read N times

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="${SCRIPT_DIR}/data"
CONFIG_JSON="${SCRIPT_DIR}/config.json"
NAIL_THREADS="${NAIL_THREADS:-$(nproc 2>/dev/null || echo 16)}"
FASTA_DIR="${SCRIPT_DIR}/fasta"
RESULTS_ROOT="${SCRIPT_DIR}/results"
OUTPUT_FILE_BASE=""
APPL=""
IPRLOOKUP=0
GOTERMS=0
BATCH_SIZE=0   # 0 = unlimited (all genomes in one call)

usage() {
  echo "Usage: $0 -b <output_file_base> [-f fasta_dir] [-o results_dir] [-t N] [--batch-size N] [--applications db1,db2] [--iprlookup] [--goterms]"
  echo "  -b            Output file base name (required); results in <results_dir>/<output_file_base>/<genome_id>.tsv.gz"
  echo "  -f            Directory containing .fa/.faa/.fasta files (default: fasta/)"
  echo "  -o            Output root directory (default: results/)"
  echo "  -t            Threads (default: nproc or 16)"
  echo "  --batch-size  Genomes per nail call (default: 0 = all).  Lower values use less RAM"
  echo "                but read the query HMM database more times.  Use when RAM is limited."
  echo "  --applications  Comma-separated list of databases to run (default: all in config)"
  echo "                  Available: pfam,ncbifam,superfamily,antifam,hamap,sfld,pirsf,pirsr,panther,cath"
  echo "  --iprlookup   Add InterPro annotation column"
  echo "  --goterms     Add GO terms column"
  echo "  --help        This message"
}
[[ " $* " == *" --help "* ]] || [[ " $* " == *" -help "* ]] && { usage; exit 0; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--output-file-base) OUTPUT_FILE_BASE="$2"; shift 2 ;;
    -f) FASTA_DIR="$2"; shift 2 ;;
    -o) RESULTS_ROOT="$2"; shift 2 ;;
    -t) NAIL_THREADS="$2"; shift 2 ;;
    --batch-size) BATCH_SIZE="$2"; shift 2 ;;
    --applications|--appl) APPL="${2:-}"; shift 2 ;;
    --iprlookup) IPRLOOKUP=1; shift ;;
    --goterms) GOTERMS=1; shift ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$OUTPUT_FILE_BASE" ]]; then
  echo "Error: -b <output_file_base> is required"
  usage
  exit 1
fi
if [[ ! -d "$FASTA_DIR" ]]; then
  echo "Error: Fasta directory not found: $FASTA_DIR"
  exit 1
fi
if ! command -v nail &>/dev/null; then
  echo "Error: nail not in PATH. Activate env: conda activate nailscan"
  exit 1
fi
if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Error: config.json not found"
  exit 1
fi

FASTA_DIR="$(cd "$FASTA_DIR" && pwd)"
RESULTS_ROOT="$(cd "$RESULTS_ROOT" 2>/dev/null || (mkdir -p "$RESULTS_ROOT" && cd "$RESULTS_ROOT" && pwd))"

DB_LIST=""
if [[ -n "$APPL" ]]; then
  DB_LIST="$(echo "$APPL" | tr ',' ' ')"
else
  DB_LIST="$(python3 -c '
import json, sys
with open(sys.argv[1]) as f:
    c = json.load(f)
print(" ".join(k for k in c if isinstance(c.get(k), dict)))
' "$CONFIG_JSON" 2>/dev/null)"
fi

# Collect ordered list of FASTA files (preserving discovery order for batch chunking)
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

FA_LIST=()
seen=""
for fa in "$FASTA_DIR"/*.fa "$FASTA_DIR"/*.faa "$FASTA_DIR"/*.fasta; do
  [[ -f "$fa" ]] || continue
  [[ " $seen " == *" $fa "* ]] && continue
  seen="$seen $fa"
  FA_LIST+=("$fa")
done
n_genomes=${#FA_LIST[@]}
if [[ $n_genomes -eq 0 ]]; then
  echo "No FASTA files found in ${FASTA_DIR}"
  exit 1
fi

# Determine chunk size: 0 (unlimited) means process all genomes in a single nail call
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

HEADER_MAP="target	NAME	ACC	DESC	target_start	target_end	query_start	query_end	score	bias	evalue	cell_frac"
HEADER_NO_MAP="target	NAME	target_start	target_end	query_start	query_end	score	bias	evalue	cell_frac"

# _run_db_on_batch <db> <hmm_path> <batch_fasta>
# Runs nail + splits results into per-genome TSVs under ${WORK}/per_genome_<db>/
_run_db_on_batch() {
  local db="$1" hmm_path="$2" batch_fasta="$3"
  local map_file="${hmm_path}.map"

  (cd "$WORK" && nail search -t "$NAIL_THREADS" --tbl-out "results_${db}.tbl" "$hmm_path" "$batch_fasta")
  if [[ ! -f "$WORK/results_${db}.tbl" ]]; then
    echo "Error: nail did not produce results for $db"
    return 1
  fi

  if [[ ! -f "$map_file" ]]; then
    awk 'BEGIN{name="";acc="";desc=""}
         /^NAME[ \t]+/ { if (name != "") print name "\t" acc "\t" desc; name=$2; acc=""; desc="" }
         /^ACC[ \t]+/  { acc=$2 }
         /^DESC[ \t]+/ { desc=$0; sub(/^DESC[ \t]+/, "", desc) }
         END { if (name != "") print name "\t" acc "\t" desc }' "$hmm_path" > "$map_file"
  fi

  mkdir -p "${WORK}/per_genome_${db}"
  if [[ -f "$map_file" ]]; then
    awk -v work="${WORK}/per_genome_${db}" -v hdr="$HEADER_MAP" -v mapf="$map_file" '
      BEGIN {
        # Read map file with explicit tab splitting to preserve multi-word descriptions
        while ((getline line < mapf) > 0) {
          t1 = index(line, "\t")
          if (t1 == 0) continue
          key = substr(line, 1, t1-1)
          map[key] = substr(line, t1+1)   # "ACC\tFULL DESC"
        }
        close(mapf)
      }
      /^#/ || NF<10 { next }
      {
        target=$1; query=$2
        idx=index(target,"|"); if (idx==0) next
        genome_id=substr(target,1,idx-1)
        gene_id=substr(target,idx+1)
        f=work "/" genome_id ".tsv"
        if (!(genome_id in seen)) { print hdr > f; seen[genome_id]=1 }
        rest=map[query]
        if (rest == "") { acc=""; desc="" } else {
          t = index(rest, "\t"); acc=substr(rest,1,t-1); desc=substr(rest,t+1) }
        print gene_id"\t"query"\t"acc"\t"desc"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 >> f
      }
    ' "$WORK/results_${db}.tbl"
  else
    awk -v work="${WORK}/per_genome_${db}" -v hdr="$HEADER_NO_MAP" '
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
    ' "$WORK/results_${db}.tbl"
  fi

  # Run DB-specific post-processing (ID cleanup, threshold filtering, best-hit).
  # Always applied for dbs with database-specific logic; also for --iprlookup/--goterms.
  if [[ $IPRLOOKUP -eq 1 ]] || [[ $GOTERMS -eq 1 ]] || \
     [[ "$db" == "panther" ]] || [[ "$db" == "hamap" ]] || \
     [[ "$db" == "superfamily" ]] || [[ "$db" == "pfam" ]] || \
     [[ "$db" == "ncbifam" ]]; then
    if [[ -f "${SCRIPT_DIR}/scripts/add_ipr_go.py" ]]; then
      for tsv in "${WORK}/per_genome_${db}"/*.tsv; do
        [[ -f "$tsv" ]] || continue
        python3 "${SCRIPT_DIR}/scripts/add_ipr_go.py" --tsv "$tsv" --db "$db" --data-dir "$DATA_ROOT" \
          $([[ $IPRLOOKUP -eq 1 ]] && echo --iprlookup) $([[ $GOTERMS -eq 1 ]] && echo --goterms) > "${tsv}.tmp"
        mv "${tsv}.tmp" "$tsv"
      done
    fi
  fi
  rm -f "$WORK/results_${db}.tbl"
}

# Build DB info once (hmm_name + path) before the chunk loop
declare -A DB_HMM_PATH
for db in $DB_LIST; do
  hmm_name="$(python3 -c '
import json, sys
with open(sys.argv[1]) as f:
    c = json.load(f)
v = c.get(sys.argv[2], {})
print(v.get("hmm_name", "") if isinstance(v, dict) else "")
' "$CONFIG_JSON" "$db" 2>/dev/null)" || true
  if [[ -z "$hmm_name" ]]; then echo "Skipping $db - no hmm_name"; continue; fi
  hmm_path="${DATA_ROOT}/${db}/${hmm_name}"
  if [[ ! -f "$hmm_path" ]]; then echo "Skipping $db - HMM not found: $hmm_path"; continue; fi
  DB_HMM_PATH["$db"]="$(cd "$(dirname "$hmm_path")" && pwd)/$(basename "$hmm_path")"
done

# Process genomes in chunks
chunk_idx=0
start=0
while [[ $start -lt $n_genomes ]]; do
  end=$(( start + chunk ))
  [[ $end -gt $n_genomes ]] && end=$n_genomes
  chunk_idx=$(( chunk_idx + 1 ))

  # Build concatenated FASTA for this chunk
  BATCH_FASTA="${WORK}/batch_${chunk_idx}.fa"
  : > "$BATCH_FASTA"
  for (( i=start; i<end; i++ )); do
    fa="${FA_LIST[$i]}"
    base="$(basename "$fa")"
    genome_id="$base"
    [[ "$base" == *.fasta ]] && genome_id="${base%.fasta}"
    [[ "$base" == *.faa   ]] && genome_id="${base%.faa}"
    [[ "$base" == *.fa    ]] && genome_id="${base%.fa}"
    awk -v g="$genome_id" '/^>/ { print ">" g "|" substr($0,2); next } { print }' "$fa" >> "$BATCH_FASTA"
  done
  n_seqs=$(grep -c '^>' "$BATCH_FASTA" || true)
  [[ $n_chunks -gt 1 ]] && echo "" && echo "=== Batch ${chunk_idx}/${n_chunks}: genomes $((start+1))-${end} (${n_seqs} sequences) ==="

  for db in "${!DB_HMM_PATH[@]}"; do
    echo ""
    echo "Starting $db ..."
    _run_db_on_batch "$db" "${DB_HMM_PATH[$db]}" "$BATCH_FASTA"
  done

  rm -f "$BATCH_FASTA"
  start=$end
done

# Merge per-DB temp TSVs into one TSV per genome under OUTPUT_FILE_BASE
OUTPUT_DIR="${RESULTS_ROOT}/${OUTPUT_FILE_BASE}"
mkdir -p "$OUTPUT_DIR"
genome_ids=""
for d in "${WORK}"/per_genome_*; do
  [[ -d "$d" ]] || continue
  for tsv in "$d"/*.tsv; do
    [[ -f "$tsv" ]] || continue
    genome_ids="$genome_ids $(basename "$tsv" .tsv)"
  done
  break
done
genome_ids="$(echo "$genome_ids" | tr ' ' '\n' | sort -u)"

for genome_id in $genome_ids; do
  [[ -z "$genome_id" ]] && continue
  first=1
  merged="${WORK}/merged_${genome_id}.tsv"
  : > "$merged"
  for d in "${WORK}"/per_genome_*; do
    [[ -d "$d" ]] || continue
    tsv="${d}/${genome_id}.tsv"
    [[ -f "$tsv" ]] || continue
    db="$(basename "$d" | sed 's/^per_genome_//')"
    if [[ $first -eq 1 ]]; then
      awk -v db="$db" 'NR==1 { print "Analysis\t" $0; next } { print db "\t" $0 }' "$tsv" >> "$merged"
      first=0
    else
      awk -v db="$db" 'NR>1 { print db "\t" $0 }' "$tsv" >> "$merged"
    fi
  done
  if [[ -s "$merged" ]]; then
    gzip -c "$merged" > "${OUTPUT_DIR}/${genome_id}.tsv.gz"
  fi
done

echo ""
echo "Done. Results: ${OUTPUT_DIR}/"
