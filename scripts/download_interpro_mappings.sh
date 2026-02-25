#!/usr/bin/env bash
# Download InterPro mapping files for -iprlookup and -goterms.
# Puts files in data/ (interpro2go, signature2interpro.tsv) and builds signature2interpro from interpro.xml.gz.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="${SCRIPT_DIR}/../data"
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release"

mkdir -p "$DATA_ROOT"
cd "$DATA_ROOT"

echo "==> Downloading interpro2go ..."
if [[ ! -f interpro2go ]] || [[ "${1:-}" == "--force" ]]; then
  curl -sL -o interpro2go "${BASE_URL}/interpro2go"
  echo "    $(wc -l < interpro2go) lines"
else
  echo "    already present (use --force to re-download)"
fi

echo "==> Downloading interpro.xml.gz (for signature->InterPro mapping) ..."
if [[ ! -f interpro.xml.gz ]] || [[ "${1:-}" == "--force" ]]; then
  curl -sL -o interpro.xml.gz "${BASE_URL}/interpro.xml.gz"
  echo "    $(ls -lh interpro.xml.gz | awk '{print $5}')"
else
  echo "    already present (use --force to re-download)"
fi

echo "==> Building signature2interpro.tsv from interpro.xml.gz ..."
python3 "${SCRIPT_DIR}/build_signature2interpro.py" interpro.xml.gz > signature2interpro.tsv.tmp
mv signature2interpro.tsv.tmp signature2interpro.tsv
echo "    $(wc -l < signature2interpro.tsv) mappings"

echo "==> InterPro/GO mapping data in ${DATA_ROOT}/"
echo "    interpro2go, signature2interpro.tsv"
