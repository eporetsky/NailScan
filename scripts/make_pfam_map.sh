#!/usr/bin/env bash
# Parse a Pfam HMM file; write NAME, ACC, DESC per model to a tab-separated map file.
# Usage: make_pfam_map.sh <pfam.hmm> [output.map]
# Default output: <pfam.hmm>.map

set -euo pipefail
HMM="${1:-}"
OUT="${2:-}"

if [[ -z "$HMM" ]] || [[ ! -f "$HMM" ]]; then
  echo "Usage: $0 <pfam.hmm> [output.map]" >&2
  exit 1
fi
[[ -n "$OUT" ]] || OUT="${HMM}.map"

awk '/^NAME[ \t]+/ { if (name != "") print name "\t" acc "\t" desc; name=$2; acc=""; desc="" }
     /^ACC[ \t]+/  { acc=$2 }
     /^DESC[ \t]+/ { desc=$0; sub(/^DESC[ \t]+/, "", desc) }
     END { if (name != "") print name "\t" acc "\t" desc }' "$HMM" > "$OUT"
echo "Wrote ${OUT} ($(wc -l < "$OUT") models)"
