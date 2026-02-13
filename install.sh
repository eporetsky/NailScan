#!/usr/bin/env bash
set -euo pipefail

# NailScan install: build nail and optionally download Pfam DB to data/pfam (from config.json)
# Run with the nail conda env activated: conda activate nail && ./install.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAIL_REPO_URL="${NAIL_REPO_URL:-https://github.com/TravisWheelerLab/nail}"
NAIL_DIR="${NAIL_DIR:-${SCRIPT_DIR}/nail-repo}"
NAIL_VERSION="${NAIL_VERSION:-main}"
DATA_PFAM="${SCRIPT_DIR}/data/pfam"
CONFIG_JSON="${SCRIPT_DIR}/config.json"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "Error: Conda environment not active. Run: conda activate nail"
  exit 1
fi

if ! command -v cargo &>/dev/null; then
  echo "Error: cargo not found. Install the nail env: conda env create -f environment.yml && conda activate nail"
  exit 1
fi

if ! command -v mmseqs &>/dev/null; then
  echo "Error: mmseqs not found in PATH. The nail conda env should provide mmseqs2."
  exit 1
fi

echo "==> Cloning nail (version: ${NAIL_VERSION}) into ${NAIL_DIR} ..."
if [[ -d "$NAIL_DIR" ]]; then
  (cd "$NAIL_DIR" && git fetch --all && git checkout "$NAIL_VERSION" && git pull --ff-only 2>/dev/null || true)
else
  git clone --depth 1 --branch "$NAIL_VERSION" "$NAIL_REPO_URL" "$NAIL_DIR"
fi

echo "==> Building nail (release) ..."
(cd "$NAIL_DIR" && cargo build --release)

NAIL_BIN="${NAIL_DIR}/target/release/nail"
if [[ ! -f "$NAIL_BIN" ]]; then
  echo "Error: Build failed, binary not found at $NAIL_BIN"
  exit 1
fi

mkdir -p "$CONDA_PREFIX/bin"
cp "$NAIL_BIN" "$CONDA_PREFIX/bin/nail"
echo "==> Installed nail to $CONDA_PREFIX/bin/nail"

nail -h >/dev/null && echo "==> Nail OK. Run 'nail search -h' for usage."

# --- Pfam database: if data/pfam/pfam_a.hmm exists, skip; else read config.json, download, md5, unpack
if [[ -f "${DATA_PFAM}/pfam_a.hmm" ]]; then
  echo "==> Pfam database already present at ${DATA_PFAM}/ (skipping download)"
else
  if [[ ! -f "$CONFIG_JSON" ]]; then
    echo "==> No config.json found; skipping Pfam download. Put Pfam HMM in data/pfam/ or add config.json."
  else
    # Parse config.json for pfam url and md5_url (minimal grep/sed, no extra deps)
    PFAM_URL="$(sed -n 's/.*"url"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$CONFIG_JSON" | head -1)"
    MD5_URL="$(sed -n 's/.*"md5_url"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$CONFIG_JSON" | head -1)"
    if [[ -z "$PFAM_URL" ]]; then
      echo "Error: config.json must contain pfam.url"
      exit 1
    fi
    [[ -n "$MD5_URL" ]] || MD5_URL="${PFAM_URL}.md5"
    mkdir -p "$DATA_PFAM"
    TAR_TMP="$(mktemp -p "$(dirname "$DATA_PFAM")" pfam.XXXXXXXX.tar.gz)"
    trap 'rm -f "$TAR_TMP"' EXIT
    echo "==> Downloading Pfam from $PFAM_URL ..."
    if command -v wget &>/dev/null; then
      wget -q -O "$TAR_TMP" "$PFAM_URL"
    elif command -v curl &>/dev/null; then
      curl -sL -o "$TAR_TMP" "$PFAM_URL"
    else
      echo "Error: need wget or curl to download Pfam"
      exit 1
    fi
    if [[ -n "$MD5_URL" ]]; then
      MD5FILE="$(mktemp)"
      if command -v wget &>/dev/null; then wget -q -O "$MD5FILE" "$MD5_URL"; else curl -sL -o "$MD5FILE" "$MD5_URL"; fi
      # EBI md5 file may be "hash  filename"; ensure we check our tarball (rewrite to expected name for md5sum -c)
      EXPECTED_HASH="$(awk '{print $1}' "$MD5FILE")"
      ACTUAL_HASH="$(md5sum < "$TAR_TMP" | awk '{print $1}')"
      rm -f "$MD5FILE"
      if [[ "$EXPECTED_HASH" != "$ACTUAL_HASH" ]]; then
        echo "Error: md5 checksum failed for Pfam tarball (expected $EXPECTED_HASH, got $ACTUAL_HASH)"
        exit 1
      fi
      echo "==> Md5 checksum OK"
    fi
    echo "==> Unpacking Pfam to ${DATA_PFAM}/ ..."
    # Tarball has pfam/38.1/pfam_a.hmm and pfam_a.dat; extract then move into data/pfam
    UNPACK_DIR="$(mktemp -d -p "$(dirname "$DATA_PFAM")" pfam_unpack.XXXXXXXX)"
    tar -xzf "$TAR_TMP" -C "$UNPACK_DIR"
    if [[ -d "$UNPACK_DIR/pfam" ]]; then
      SUB="$(find "$UNPACK_DIR/pfam" -mindepth 1 -maxdepth 1 -type d | head -1)"
      if [[ -n "$SUB" ]] && [[ -f "${SUB}/pfam_a.hmm" ]]; then
        mv "${SUB}"/* "$DATA_PFAM/"
      else
        find "$UNPACK_DIR/pfam" -type f -exec mv {} "$DATA_PFAM/" \;
      fi
    fi
    rm -rf "$UNPACK_DIR"
    rm -f "$TAR_TMP"
    trap - EXIT
    if [[ -f "${DATA_PFAM}/pfam_a.hmm" ]]; then
      echo "==> Pfam database installed at ${DATA_PFAM}/"
    else
      echo "Error: Pfam unpack did not produce data/pfam/pfam_a.hmm"
      exit 1
    fi
  fi
fi

# Generate Pfam NAME -> ACC,DESC mapping. After install: source install.sh && make_pfam_map data/pfam/pfam_a.hmm
make_pfam_map() {
  local hmm="${1:-data/pfam/pfam_a.hmm}"
  local out="${2:-${hmm}.map}"
  [[ -f "$hmm" ]] || { echo "HMM not found: $hmm"; return 1; }
  awk '/^NAME[ \t]+/ { if (name != "") print name "\t" acc "\t" desc; name=$2; acc=""; desc="" }
       /^ACC[ \t]+/  { acc=$2 }
       /^DESC[ \t]+/ { desc=$0; sub(/^DESC[ \t]+/, "", desc) }
       END { if (name != "") print name "\t" acc "\t" desc }' "$hmm" > "$out"
  echo "Wrote ${out} ($(wc -l < "$out") models)"
}
