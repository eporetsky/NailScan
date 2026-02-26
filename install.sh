#!/usr/bin/env bash
set -euo pipefail

# NailScan install: build nail and optionally download HMM DBs to data/<dbname> (from config.json).
# Run with the nailscan conda env activated: conda activate nailscan && ./install.sh
# To install only specific DBs: NAILSCAN_DBS="pfam ncbifam" ./install.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAIL_REPO_URL="${NAIL_REPO_URL:-https://github.com/TravisWheelerLab/nail}"
NAIL_DIR="${NAIL_DIR:-${SCRIPT_DIR}/nail-repo}"
# Pinned to the commit the binarize patch was written against.
# To upgrade nail: update this hash and regenerate patches/nail-binarize.patch.
NAIL_VERSION="${NAIL_VERSION:-b6291d398d4af4dddc06fd8783b0b95990ca8e49}"
DATA_ROOT="${SCRIPT_DIR}/data"
CONFIG_JSON="${SCRIPT_DIR}/config.json"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "Error: Conda environment not active. Run: conda activate nailscan"
  exit 1
fi

if ! command -v cargo &>/dev/null; then
  echo "Error: cargo not found. Install the nailscan env: conda env create -f environment.yml && conda activate nailscan"
  exit 1
fi

if ! command -v mmseqs &>/dev/null; then
  echo "Error: mmseqs not found in PATH. The nail conda env should provide mmseqs2."
  exit 1
fi

# --- Build nail ---
echo "==> Cloning nail (version: ${NAIL_VERSION}) into ${NAIL_DIR} ..."
if [[ -d "$NAIL_DIR" ]]; then
  (cd "$NAIL_DIR" && git fetch --all && git checkout "$NAIL_VERSION" && git pull --ff-only 2>/dev/null || true)
else
  git clone --depth 1 --branch "$NAIL_VERSION" "$NAIL_REPO_URL" "$NAIL_DIR"
fi

PATCH_FILE="${SCRIPT_DIR}/patches/nail-binarize.patch"
if [[ -f "$PATCH_FILE" ]]; then
  echo "==> Applying nailscan patches ..."
  (cd "$NAIL_DIR" && git apply "$PATCH_FILE")
else
  echo "Warning: patches/nail-binarize.patch not found; building unpatched nail (DB loading will be slower)"
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

if command -v nail &>/dev/null; then
  NAIL_WHICH="$(which nail)"
  if [[ "$NAIL_WHICH" != "$CONDA_PREFIX/bin/nail" ]]; then
    echo "Warning: 'nail' in PATH is $NAIL_WHICH (not the one just installed). Run: $CONDA_PREFIX/bin/nail"
  fi
fi
nail -h >/dev/null && echo "==> Nail OK. Run 'nail search -h' for usage."

# --- Install HMM databases from config.json ---
# One script for all DBs: download + MD5 + unpack are generic; placement into data/<dbname>/
# is done by place_db_files(). Add a case there for any DB that needs different processing
# (e.g. multiple HMM files, conversion step). Default: find file named hmm_name, move its
# containing directory's contents into data_dir.
place_db_files() {
  local dbname="$1"
  local unpack_dir="$2"
  local data_dir="$3"
  local hmm_name="$4"
  case "$dbname" in
    # Add DB-specific placement here, e.g.:
    # mydb)  cp -r "$unpack_dir"/special/layout/* "$data_dir"/ ;;
    *)
      # Default: one HMM file (or one concatenated file); find it and move its directory contents
      local hmm_path
      hmm_path="$(find "$unpack_dir" -name "$hmm_name" -type f | head -1)"
      if [[ -z "$hmm_path" ]]; then
        echo "Error: ${dbname} tarball did not contain a file named '${hmm_name}'"
        return 1
      fi
      local hmm_dir
      hmm_dir="$(dirname "$hmm_path")"
      mv "$hmm_dir"/* "$data_dir/" 2>/dev/null || true
      ;;
  esac
}

install_one_db() {
  local dbname="$1"
  local url="$2"
  local md5_url="$3"
  local hmm_name="$4"
  local data_dir="${DATA_ROOT}/${dbname}"

  if [[ -f "${data_dir}/${hmm_name}" ]]; then
    echo "==> ${dbname} database already present at ${data_dir}/ (skipping download)"
    return 0
  fi

  if [[ -z "$url" ]]; then
    echo "Error: config.json must contain ${dbname}.url"
    return 1
  fi

  [[ -n "$md5_url" ]] || md5_url="${url}.md5"
  mkdir -p "$data_dir"
  local tar_tmp
  tar_tmp="$(mktemp -p "$(dirname "$data_dir")" "${dbname}.XXXXXXXX.tar.gz")"
  trap 'rm -f "$tar_tmp"' EXIT

  echo "==> Downloading ${dbname} from ${url} ..."
  if command -v wget &>/dev/null; then
    wget -q -O "$tar_tmp" "$url"
  elif command -v curl &>/dev/null; then
    curl -sL -o "$tar_tmp" "$url"
  else
    echo "Error: need wget or curl to download ${dbname}"
    return 1
  fi

  if [[ -n "$md5_url" ]]; then
    local md5file
    md5file="$(mktemp)"
    if command -v wget &>/dev/null; then wget -q -O "$md5file" "$md5_url"; else curl -sL -o "$md5file" "$md5_url"; fi
    local expected_hash actual_hash
    expected_hash="$(awk '{print $1}' "$md5file")"
    actual_hash="$(md5sum < "$tar_tmp" | awk '{print $1}')"
    rm -f "$md5file"
    if [[ "$expected_hash" != "$actual_hash" ]]; then
      echo "Error: md5 checksum failed for ${dbname} tarball (expected $expected_hash, got $actual_hash)"
      return 1
    fi
    echo "==> Md5 checksum OK (${dbname})"
  fi

  echo "==> Unpacking ${dbname} to ${data_dir}/ ..."
  local unpack_dir
  unpack_dir="$(mktemp -d -p "$(dirname "$data_dir")" "${dbname}_unpack.XXXXXXXX")"
  tar -xzf "$tar_tmp" -C "$unpack_dir"

  if ! place_db_files "$dbname" "$unpack_dir" "$data_dir" "$hmm_name"; then
    rm -rf "$unpack_dir"
    rm -f "$tar_tmp"
    trap - EXIT
    return 1
  fi

  rm -rf "$unpack_dir"
  rm -f "$tar_tmp"
  trap - EXIT

  if [[ -f "${data_dir}/${hmm_name}" ]]; then
    echo "==> ${dbname} database installed at ${data_dir}/"
  else
    echo "Error: ${dbname} unpack did not produce ${data_dir}/${hmm_name}"
    return 1
  fi
}

if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "==> No config.json found; skipping database downloads. Add config.json to install DBs into data/."
else
  # Optional: install only these DBs (default: all in config)
  want_dbs="${NAILSCAN_DBS:-}"
  while IFS= read -r line; do
    [[ -z "$line" ]] && continue
    dbname="${line%%	*}"
    rest="${line#*	}"
    url="${rest%%	*}"
    rest="${rest#*	}"
    md5_url="${rest%%	*}"
    hmm_name="${rest#*	}"
    if [[ -n "$want_dbs" ]] && [[ " $want_dbs " != *" $dbname "* ]]; then
      continue
    fi
    install_one_db "$dbname" "$url" "$md5_url" "$hmm_name"
  done < <(python3 -c '
import json, sys
with open(sys.argv[1]) as f:
    c = json.load(f)
for k, v in c.items():
    if not isinstance(v, dict):
        continue
    url = v.get("url", "")
    md5_url = v.get("md5_url", "")
    hmm_name = v.get("hmm_name", "")
    if not hmm_name and url:
        if "pfam" in k.lower():
            hmm_name = "pfam_a.hmm"
        elif "ncbifam" in k.lower():
            hmm_name = "ncbifam.hmm"
        elif "superfamily" in k.lower():
            hmm_name = "hmmlib_1.75"
        elif "antifam" in k.lower():
            hmm_name = "AntiFam.hmm"
        elif "hamap" in k.lower():
            hmm_name = "hamap.hmm.lib"
        elif "sfld" in k.lower():
            hmm_name = "sfld.hmm"
        elif "pirsf" in k.lower():
            hmm_name = "sf_hmm_all"
        elif "pirsr" in k.lower():
            hmm_name = "sr_hmm_all"
        elif "panther" in k.lower():
            hmm_name = "binHmm"
        elif "cath" in k.lower():
            hmm_name = "gene3d_main.hmm"
    print(k, url, md5_url, hmm_name, sep="\t")
' "$CONFIG_JSON")
fi

# --- Fix HMMER2 / HMMER3/b databases: convert to HMMER3/f + add COMPO lines ---
# nail's from_p7hmm parser requires a COMPO line in each model's header to advance
# its state machine into the model body.  COMPO is optional per the HMMER3 spec but
# nail requires it.  Some databases (Superfamily, PIRSR) are distributed in HMMER2
# binary format (reported as "HMMER3/b"); hmmconvert produces valid HMMER3/f but
# omits COMPO lines, so all models parse as empty and nail returns 0 hits.
# This function converts + inserts COMPO lines for any HMM file that needs it.
# The original file is saved as <hmm>.raw so the fix is skipped on re-runs.
# See docs/NAIL_PATCH.md for the full explanation.
_fix_hmm3b_if_needed() {
  local hmm_file="$1"
  local hmm_backup="${hmm_file}.raw"
  local db_dir
  db_dir="$(dirname "$hmm_file")"

  [[ -f "$hmm_file" ]] || return 0

  # Skip if already processed (backup exists AND file has COMPO lines)
  if [[ -f "$hmm_backup" ]] && grep -qm1 "^  COMPO" "$hmm_file" 2>/dev/null; then
    echo "==> $(basename "$hmm_file"): already nail-ready (HMMER3+COMPO); skipping"
    return 0
  fi

  local first_line
  first_line="$(head -1 "$hmm_file" 2>/dev/null || true)"
  local needs_convert=false
  if [[ "$first_line" == *"HMMER3/b"* ]] || [[ "$first_line" == *"HMMER2"* ]]; then
    needs_convert=true
  fi

  # If already HMMER3/f and has COMPO, nothing to do (e.g. no backup yet but fine)
  if [[ "$needs_convert" == false ]] && grep -qm1 "^  COMPO" "$hmm_file" 2>/dev/null; then
    return 0
  fi

  # Determine the source to convert from
  local source
  if [[ -f "$hmm_backup" ]]; then
    echo "==> $(basename "$hmm_file"): re-processing from backup (adding COMPO lines) ..."
    source="$hmm_backup"
  else
    echo "==> $(basename "$hmm_file"): processing for nail (HMMER2→HMMER3 + COMPO lines) ..."
    cp "$hmm_file" "$hmm_backup"
    source="$hmm_backup"
  fi

  local conv_tmp
  conv_tmp="$(mktemp -p "$db_dir" hmm_conv.XXXXXXXX)"

  if [[ "$needs_convert" == true ]]; then
    if ! command -v hmmconvert &>/dev/null; then
      echo "Warning: hmmconvert not found (conda install hmmer); skipping $(basename "$hmm_file") fix."
      rm -f "$conv_tmp"
      return 1
    fi
    if ! hmmconvert "$source" > "$conv_tmp" 2>/dev/null; then
      echo "Warning: hmmconvert failed for $(basename "$hmm_file")."
      rm -f "$conv_tmp"
      return 1
    fi
  else
    cp "$source" "$conv_tmp"
  fi

  # Insert a COMPO line after each transition header ("m->m ...") when one is absent.
  # Values are taken from the position-0 insert-emission line already present, which
  # holds the standard HMMER3 null-model amino-acid background frequencies.
  local final_tmp
  final_tmp="$(mktemp -p "$db_dir" hmm_nail.XXXXXXXX)"
  python3 - "$conv_tmp" "$final_tmp" <<'PYEOF'
import sys
infile, outfile = sys.argv[1], sys.argv[2]
with open(infile) as f:
    lines = f.readlines()
result = []
i, n = 0, len(lines)
while i < n:
    line = lines[i]
    tokens = line.split()
    if tokens and tokens[0] == 'm->m':
        result.append(line)
        i += 1
        while i < n and not lines[i].strip():
            result.append(lines[i])
            i += 1
        if i < n:
            nt = lines[i].split()
            if nt and nt[0] != 'COMPO' and len(nt) >= 20:
                result.append('  COMPO   ' + '  '.join(nt[:20]) + '\n')
    else:
        result.append(line)
        i += 1
with open(outfile, 'w') as f:
    f.writelines(result)
PYEOF

  rm -f "$conv_tmp"
  mv "$final_tmp" "$hmm_file"
  echo "==> $(basename "$hmm_file"): nail-ready HMM installed (original → $(basename "$hmm_backup"))"
}

# Apply to all installed databases (no-op if already processed or not HMMER2)
if [[ -f "$CONFIG_JSON" ]]; then
  while IFS= read -r line; do
    [[ -z "$line" ]] && continue
    dbname="${line%%	*}"
    hmm_name="${line##*	}"
    hmm_file="${DATA_ROOT}/${dbname}/${hmm_name}"
    [[ -f "$hmm_file" ]] && _fix_hmm3b_if_needed "$hmm_file"
  done < <(python3 -c '
import json, sys
with open(sys.argv[1]) as f:
    c = json.load(f)
for k, v in c.items():
    if isinstance(v, dict) and v.get("hmm_name"):
        print(k, v["hmm_name"], sep="\t")
' "$CONFIG_JSON")
fi

# --- Pre-binarize HMM databases for fast loading ---
# nail prep parses each HMM file and writes a <hmm>.nailbin bincode file that
# nail search loads in ~2s instead of ~56s for large databases like Pfam.
if [[ -f "$CONFIG_JSON" ]] && command -v nail &>/dev/null; then
  while IFS= read -r line; do
    [[ -z "$line" ]] && continue
    dbname="${line%%	*}"
    hmm_name="${line##*	}"
    hmm_file="${DATA_ROOT}/${dbname}/${hmm_name}"
    nailbin_file="${hmm_file}.nailbin"
    if [[ ! -f "$hmm_file" ]]; then
      continue
    fi
    if [[ -f "$nailbin_file" ]] && [[ "$nailbin_file" -nt "$hmm_file" ]]; then
      echo "==> ${dbname}: binary DB already up to date (skipping prep)"
      continue
    fi
    echo "==> Binarizing ${dbname} query database ..."
    nail prep "$hmm_file"
  done < <(python3 -c '
import json, sys
with open(sys.argv[1]) as f:
    c = json.load(f)
for k, v in c.items():
    if isinstance(v, dict) and v.get("hmm_name"):
        print(k, v["hmm_name"], sep="\t")
' "$CONFIG_JSON")
fi

# --- PANTHER family names (for DESC column and post-processing) ---
# Downloaded from data.pantherdb.org; used by scripts/add_ipr_go.py to add human-readable
# family descriptions and to apply the best-hit filter that InterProScan uses.
PANTHER_DIR="${DATA_ROOT}/panther"
if [[ -d "$PANTHER_DIR" ]] && ! ls "${PANTHER_DIR}"/PANTHER*_HMM_classifications 2>/dev/null | grep -q .; then
  PANTHER_NAMES_URL="https://data.pantherdb.org/ftp/hmm_classifications/current_release/PANTHER19.0_HMM_classifications"
  echo "==> Downloading PANTHER family names ..."
  if command -v wget &>/dev/null; then
    wget -q -O "${PANTHER_DIR}/PANTHER19.0_HMM_classifications" "$PANTHER_NAMES_URL" && \
      echo "==> PANTHER family names downloaded" || \
      echo "Warning: PANTHER names download failed (DESC column will be empty)"
  elif command -v curl &>/dev/null; then
    curl -sL -o "${PANTHER_DIR}/PANTHER19.0_HMM_classifications" "$PANTHER_NAMES_URL" && \
      echo "==> PANTHER family names downloaded" || \
      echo "Warning: PANTHER names download failed (DESC column will be empty)"
  fi
elif [[ -d "$PANTHER_DIR" ]]; then
  echo "==> PANTHER family names already present"
fi

# --- InterPro / GO mapping (for --iprlookup and --goterms) ---
SIGNATURE2IPR="${DATA_ROOT}/signature2interpro.tsv"
INTERPRO2GO="${DATA_ROOT}/interpro2go"
if [[ ! -f "$SIGNATURE2IPR" ]] || [[ ! -f "$INTERPRO2GO" ]]; then
  echo "==> InterPro/GO mapping missing; running scripts/download_interpro_mappings.sh ..."
  "${SCRIPT_DIR}/scripts/download_interpro_mappings.sh"
else
  echo "==> InterPro/GO mapping already present (signature2interpro.tsv, interpro2go)"
fi

# Generate NAME -> ACC,DESC mapping for any HMM. After install: source install.sh && make_pfam_map data/pfam/pfam_a.hmm
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

# --- Threshold lookup tables (GA for Pfam, TC for NCBIfam) ---
# Pre-build small tab-separated files used by scripts/add_ipr_go.py for filtering.
# Format: NAME<tab>SEQ_THRESHOLD<tab>DOM_THRESHOLD
_build_thresholds() {
  local hmm_path="$1"
  local out_path="$2"
  local field="$3"   # GA or TC
  [[ -f "$hmm_path" ]] || return 0
  [[ -f "$out_path" ]] && return 0   # already built
  echo "==> Building $field threshold table: $out_path ..."
  python3 - "$hmm_path" "$out_path" "$field" <<'PYEOF'
import sys
hmm_path, out_path, field = sys.argv[1], sys.argv[2], sys.argv[3]
with open(hmm_path) as f, open(out_path, 'w') as out:
    name = None
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('NAME ') or line.startswith('NAME\t'):
            name = line.split(None, 1)[1].strip()
        elif name and (line.startswith(field + ' ') or line.startswith(field + '\t')):
            parts = line.split()
            if len(parts) >= 3:
                seq_t = parts[1]
                dom_t = parts[2].rstrip(';')
                out.write(f'{name}\t{seq_t}\t{dom_t}\n')
            name = None
PYEOF
  echo "==> Done: $out_path ($(wc -l < "$out_path") models)"
}

PFAM_HMM="${DATA_ROOT}/pfam/pfam_a.hmm"
PFAM_GA="${DATA_ROOT}/pfam/pfam_a.ga"
if [[ -f "$PFAM_HMM" ]]; then
  _build_thresholds "$PFAM_HMM" "$PFAM_GA" "GA"
fi

NCBIFAM_HMM="${DATA_ROOT}/ncbifam/ncbifam.hmm"
NCBIFAM_TC="${DATA_ROOT}/ncbifam/ncbifam.tc"
if [[ -f "$NCBIFAM_HMM" ]]; then
  _build_thresholds "$NCBIFAM_HMM" "$NCBIFAM_TC" "TC"
fi

# --- SFLD GA threshold table ---
SFLD_HMM="${DATA_ROOT}/sfld/sfld.hmm"
SFLD_GA="${DATA_ROOT}/sfld/sfld.ga"
if [[ -f "$SFLD_HMM" ]]; then
  _build_thresholds "$SFLD_HMM" "$SFLD_GA" "GA"
fi

# --- PIRSF score threshold table (from pirsf.dat: avg - 3.5 * std_score) ---
PIRSF_DAT="${DATA_ROOT}/pirsf/pirsf.dat"
PIRSF_THRESH="${DATA_ROOT}/pirsf/pirsf.thresh"
if [[ -f "$PIRSF_DAT" ]] && [[ ! -f "$PIRSF_THRESH" ]]; then
  echo "==> Building PIRSF threshold table: $PIRSF_THRESH ..."
  python3 - "$PIRSF_DAT" "$PIRSF_THRESH" <<'PYEOF'
import sys
dat_path, out_path = sys.argv[1], sys.argv[2]
with open(dat_path) as f:
    lines = [l.rstrip('\n') for l in f]
with open(out_path, 'w') as out:
    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            pirsf_id = lines[i][1:].strip()
            if i + 2 < len(lines):
                nums = lines[i + 2].split()
                if len(nums) >= 2:
                    try:
                        avg = float(nums[0])
                        std = float(nums[1])
                        threshold = max(0.0, avg - 3.5 * std)
                        out.write(f'{pirsf_id}\t{threshold:.4f}\t{threshold:.4f}\n')
                    except ValueError:
                        pass
            i += 4
        else:
            i += 1
PYEOF
  echo "==> Done: $PIRSF_THRESH ($(wc -l < "$PIRSF_THRESH") families)"
elif [[ -f "$PIRSF_THRESH" ]]; then
  echo "==> PIRSF threshold table already present"
fi

# --- InterPro entry name lookup table (IPR_ID -> short description) ---
# Parsed from interpro.xml.gz (downloaded by scripts/download_interpro_mappings.sh).
# Used by add_ipr_go.py to populate the IPR_desc output column.
INTERPRO_XML="${DATA_ROOT}/interpro.xml.gz"
INTERPRO_NAMES="${DATA_ROOT}/interpro_names.tsv"
if [[ -f "$INTERPRO_XML" ]] && [[ ! -f "$INTERPRO_NAMES" ]]; then
  echo "==> Building interpro_names.tsv from interpro.xml.gz ..."
  python3 - "$INTERPRO_XML" "$INTERPRO_NAMES" <<'PYEOF'
import sys, gzip, re
xml_path, out_path = sys.argv[1], sys.argv[2]
id_re = re.compile(r'<interpro\s+id="(IPR\d+)"')
name_re = re.compile(r'<name>(.*?)</name>')
opener = gzip.open if xml_path.endswith('.gz') else open
with opener(xml_path, 'rt', encoding='utf-8') as f, open(out_path, 'w') as out:
    current_id = None
    for line in f:
        if current_id is None:
            m = id_re.search(line)
            if m:
                current_id = m.group(1)
        else:
            m = name_re.search(line)
            if m:
                out.write(f'{current_id}\t{m.group(1)}\n')
                current_id = None
            elif '<interpro ' in line:
                current_id = None  # missed the name, reset
PYEOF
  echo "==> Done: $INTERPRO_NAMES ($(wc -l < "$INTERPRO_NAMES") entries)"
elif [[ -f "$INTERPRO_NAMES" ]]; then
  echo "==> interpro_names.tsv already present"
fi
