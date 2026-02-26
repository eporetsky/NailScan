#!/usr/bin/env bash
# NailScan single: run nail HMM search on one FASTA; by default run all DBs from config.
# Usage: ./nailscan.single.sh -f <input.fasta> -b <output_file_base> [-o output_dir] [-t N] [--applications ...] [--iprlookup] [--goterms]
# All DB outputs are merged into one TSV: <output_dir>/<output_file_base>.tsv

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="${SCRIPT_DIR}/data"
CONFIG_JSON="${SCRIPT_DIR}/config.json"
NAIL_THREADS="${NAIL_THREADS:-$(nproc 2>/dev/null || echo 16)}"
INPUT_FASTA=""
OUTPUT_DIR="${SCRIPT_DIR}/results"
OUTPUT_FILE_BASE=""
SINGLE_HMM=""
APPL=""
IPRLOOKUP=0
GOTERMS=0

usage() {
  echo "Usage: $0 -f <input.fasta> -b <output_file_base> [-h hmm_path] [-o output_dir] [-t N] [--applications db1,db2] [--iprlookup] [--goterms]"
  echo "  -f  FASTA file (required)"
  echo "  -b  Output file base name (required); result is <output_dir>/<output_file_base>.tsv"
  echo "  -h  Single HMM path (optional; if set, run only this HMM)"
  echo "  -o  Output directory (default: results/)"
  echo "  -t  Threads (default: nproc or 16)"
  echo "  --applications  Comma-separated list of databases to run (default: all in config)"
  echo "                  Available: pfam,ncbifam,superfamily,antifam,hamap,sfld,pirsf,pirsr,panther,cath"
  echo "  --iprlookup     Add InterPro annotation column (requires data/signature2interpro.tsv)"
  echo "  --goterms       Add GO terms column (requires data/interpro2go)"
  echo "  --help          This message"
}
[[ " $* " == *" --help "* ]] || [[ " $* " == *" -help "* ]] && { usage; exit 0; }

# Parse short and long options
while [[ $# -gt 0 ]]; do
  case "$1" in
    -f) INPUT_FASTA="$2"; shift 2 ;;
    -b|--output-file-base) OUTPUT_FILE_BASE="$2"; shift 2 ;;
    -h) SINGLE_HMM="$2"; shift 2 ;;
    -o) OUTPUT_DIR="$2"; shift 2 ;;
    -t) NAIL_THREADS="$2"; shift 2 ;;
    --applications|--appl) APPL="${2:-}"; shift 2 ;;
    --iprlookup) IPRLOOKUP=1; shift ;;
    --goterms) GOTERMS=1; shift ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$INPUT_FASTA" ]]; then
  echo "Error: -f <input.fasta> is required"
  usage
  exit 1
fi
if [[ -z "$OUTPUT_FILE_BASE" ]]; then
  echo "Error: -b <output_file_base> is required"
  usage
  exit 1
fi
if [[ ! -f "$INPUT_FASTA" ]]; then
  echo "Error: Fasta file not found: $INPUT_FASTA"
  exit 1
fi
if ! command -v nail &>/dev/null; then
  echo "Error: nail not in PATH. Activate env: conda activate nailscan"
  exit 1
fi

INPUT_FASTA="$(cd "$(dirname "$INPUT_FASTA")" && pwd)/$(basename "$INPUT_FASTA")"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
base="$(basename "$INPUT_FASTA")"
genome_id="${base%.*}"
[[ "$base" == *.fasta ]] && genome_id="${base%.fasta}"
[[ "$base" == *.faa ]] && genome_id="${base%.faa}"
[[ "$base" == *.fa ]] && genome_id="${base%.fa}"

run_one_hmm() {
  local hmm_path="$1"
  local db_label="$2"
  local out_tsv="$3"
  local work_dir="$4"
  local query_path="$hmm_path"
  [[ -f "${hmm_path}.nailbin" ]] && query_path="${hmm_path}.nailbin"
  (cd "$work_dir" && nail search -t "$NAIL_THREADS" --tbl-out results.tbl "$query_path" "$INPUT_FASTA")
  if [[ ! -f "$work_dir/results.tbl" ]]; then
    echo "Error: nail did not produce results.tbl for $db_label"
    return 1
  fi
  local map_file="${hmm_path}.map"
  if [[ ! -f "$map_file" ]]; then
    awk 'BEGIN{name="";acc="";desc=""}
         /^NAME[ \t]+/ { if (name != "") print name "\t" acc "\t" desc; name=$2; acc=""; desc="" }
         /^ACC[ \t]+/  { acc=$2 }
         /^DESC[ \t]+/ { desc=$0; sub(/^DESC[ \t]+/, "", desc) }
         END { if (name != "") print name "\t" acc "\t" desc }' "$hmm_path" > "$map_file"
  fi
  local HEADER_MAP="target	NAME	ACC	DESC	target_start	target_end	query_start	query_end	score	bias	evalue	cell_frac"
  local HEADER_NO_MAP="target	NAME	target_start	target_end	query_start	query_end	score	bias	evalue	cell_frac"
  if [[ -f "$map_file" ]]; then
    awk -v hdr="$HEADER_MAP" -v mapf="$map_file" '
      BEGIN {
        # Read map file with explicit tab splitting to preserve multi-word descriptions
        while ((getline line < mapf) > 0) {
          t1 = index(line, "\t")
          if (t1 == 0) continue
          key = substr(line, 1, t1-1)
          rest = substr(line, t1+1)
          map[key] = rest   # rest = "ACC\tFULL DESC"
        }
        close(mapf)
        print hdr
      }
      /^#/ || NF<10 { next }
      { query=$2; rest=map[query]
        if (rest == "") { acc=""; desc="" } else {
          t = index(rest, "\t"); acc=substr(rest,1,t-1); desc=substr(rest,t+1) }
        print $1"\t"query"\t"acc"\t"desc"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }
    ' "$work_dir/results.tbl" > "$out_tsv"
  else
    awk -v hdr="$HEADER_NO_MAP" '
      BEGIN { print hdr }
      /^#/ || NF<10 { next }
      { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }
    ' "$work_dir/results.tbl" > "$out_tsv"
  fi
  # Run DB-specific post-processing (ID cleanup, threshold filtering, best-hit).
  # Always applied for dbs with database-specific logic; also for --iprlookup/--goterms.
  if [[ $IPRLOOKUP -eq 1 ]] || [[ $GOTERMS -eq 1 ]] || \
     [[ "$db_label" == "panther" ]] || [[ "$db_label" == "hamap" ]] || \
     [[ "$db_label" == "superfamily" ]] || [[ "$db_label" == "pfam" ]] || \
     [[ "$db_label" == "ncbifam" ]] || [[ "$db_label" == "sfld" ]] || \
     [[ "$db_label" == "pirsf" ]] || [[ "$db_label" == "cath" ]]; then
    if [[ -n "$db_label" ]] && [[ -f "${SCRIPT_DIR}/scripts/add_ipr_go.py" ]]; then
      local tmp_tsv="${out_tsv}.tmp"
      python3 "${SCRIPT_DIR}/scripts/add_ipr_go.py" --tsv "$out_tsv" --db "$db_label" --data-dir "$DATA_ROOT" \
        $([[ $IPRLOOKUP -eq 1 ]] && echo --iprlookup) $([[ $GOTERMS -eq 1 ]] && echo --goterms) > "$tmp_tsv"
      mv "$tmp_tsv" "$out_tsv"
    fi
  fi
}

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if [[ -n "$SINGLE_HMM" ]]; then
  if [[ ! -f "$SINGLE_HMM" ]]; then
    echo "Error: HMM not found: $SINGLE_HMM"
    exit 1
  fi
  SINGLE_HMM="$(cd "$(dirname "$SINGLE_HMM")" && pwd)/$(basename "$SINGLE_HMM")"
  SINGLE_DB=""
  if [[ "$SINGLE_HMM" == *"/data/"* ]]; then
    rest="${SINGLE_HMM#*"/data/"}"
    SINGLE_DB="${rest%%/*}"
  fi
  SINGLE_DB="${SINGLE_DB:-unknown}"
  TEMP_TSV="${WORK}/single.tsv"
  START_TIME="$(date +%s)"
  echo "Running: $INPUT_FASTA (single HMM, threads: $NAIL_THREADS)"
  echo "Starting $SINGLE_DB ..."
  run_one_hmm "$SINGLE_HMM" "$SINGLE_DB" "$TEMP_TSV" "$WORK"
  FINAL_TSV="${OUTPUT_DIR}/${OUTPUT_FILE_BASE}.tsv"
  awk -v db="$SINGLE_DB" -v lmap="pfam=Pfam,cath=Gene3D,superfamily=SUPERFAMILY,panther=PANTHER,ncbifam=NCBIfam,hamap=Hamap,sfld=SFLD,pirsf=PIRSF,pirsr=PIRSR,antifam=AntiFam" '
    BEGIN { n=split(lmap,p,","); for(i=1;i<=n;i++){split(p[i],kv,"=");m[kv[1]]=kv[2]}; label=(db in m)?m[db]:db }
    NR==1 { rest=substr($0,index($0,"\t")); print $1 "\tAnalysis" rest; next }
    { rest=substr($0,index($0,"\t")); print $1 "\t" label rest }
  ' "$TEMP_TSV" > "$FINAL_TSV"
  END_TIME="$(date +%s)"
  M=$(((END_TIME - START_TIME) / 60))
  echo "completed in ${M} minutes."
  echo "$FINAL_TSV"
  exit 0
fi

# All-DBs mode: get DB list from config (skips entries with "enabled": false)
if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Error: config.json not found; cannot run all databases."
  exit 1
fi
DB_LIST=""
if [[ -n "$APPL" ]]; then
  DB_LIST="$(echo "$APPL" | tr ',' ' ')"
else
  DB_LIST="$(python3 -c '
import json, sys
with open(sys.argv[1]) as f:
    c = json.load(f)
print(" ".join(k for k in c if isinstance(c.get(k), dict) and c[k].get("enabled", True)))
' "$CONFIG_JSON" 2>/dev/null)"
fi

# Awk snippet: reorder columns (target first, Analysis second) and apply canonical DB label
_merge_awk() {
  local db="$1"
  awk -v db="$db" -v lmap="pfam=Pfam,cath=Gene3D,superfamily=SUPERFAMILY,panther=PANTHER,ncbifam=NCBIfam,hamap=Hamap,sfld=SFLD,pirsf=PIRSF,pirsr=PIRSR,antifam=AntiFam" '
    BEGIN { n=split(lmap,p,","); for(i=1;i<=n;i++){split(p[i],kv,"=");m[kv[1]]=kv[2]}; label=(db in m)?m[db]:db }
    { rest=substr($0,index($0,"\t")); print $1 "\t" label rest }
  '
}

START_TIME="$(date +%s)"
echo "Running: $INPUT_FASTA (threads: $NAIL_THREADS) for: $DB_LIST"
first=1
for db in $DB_LIST; do
  hmm_name="$(python3 -c '
import json, sys
with open(sys.argv[1]) as f:
    c = json.load(f)
v = c.get(sys.argv[2], {})
print(v.get("hmm_name", "") if isinstance(v, dict) else "")
' "$CONFIG_JSON" "$db" 2>/dev/null)" || true
  if [[ -z "$hmm_name" ]]; then echo "Skipping $db - no hmm_name in config"; continue; fi
  HMM_PATH="${DATA_ROOT}/${db}/${hmm_name}"
  if [[ ! -f "$HMM_PATH" ]]; then
    echo "Skipping $db - HMM not found: $HMM_PATH"
    continue
  fi
  WORK_DB="${WORK}/${db}"
  mkdir -p "$WORK_DB"
  TEMP_TSV="${WORK_DB}/results.tsv"
  echo "Starting $db ..."
  run_one_hmm "$HMM_PATH" "$db" "$TEMP_TSV" "$WORK_DB"
  if [[ ! -f "$TEMP_TSV" ]]; then continue; fi
  if [[ $first -eq 1 ]]; then
    # Write header (target first, then Analysis)
    head -1 "$TEMP_TSV" | awk -v lmap="pfam=Pfam" 'BEGIN{} { rest=substr($0,index($0,"\t")); print $1 "\tAnalysis" rest }' > "${WORK}/merged.tsv"
    first=0
  fi
  tail -n +2 "$TEMP_TSV" | _merge_awk "$db" >> "${WORK}/merged.tsv"
done
FINAL_TSV="${OUTPUT_DIR}/${OUTPUT_FILE_BASE}.tsv"
if [[ -f "${WORK}/merged.tsv" ]]; then
  # Sort by gene ID (column 1 = target), keeping header on top
  { head -1 "${WORK}/merged.tsv"; tail -n +2 "${WORK}/merged.tsv" | sort -t$'\t' -k1,1; } > "$FINAL_TSV"
  END_TIME="$(date +%s)"
  M=$(((END_TIME - START_TIME) / 60))
  echo "completed in ${M} minutes."
  echo "$FINAL_TSV"
else
  END_TIME="$(date +%s)"
  M=$(((END_TIME - START_TIME) / 60))
  echo "completed in ${M} minutes. No results produced."
fi
