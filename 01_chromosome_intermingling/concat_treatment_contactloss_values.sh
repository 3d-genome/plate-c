#!/bin/bash
# concat_treatment_contactloss_values.sh
#
# Usage:
#   ./concat_treatment_contactloss_values.sh aux_data/<experiment_dir>/<treatment_file>.txt
# For each old path listed in treatment_file:
# - map old_abs_path -> new_name via rename*.txt in the same directory
# - read processed/<new_name>/contacts_unisex.info
# - take 5th column (col5)
# - compute value = 100 - col5
#
# Output:
#   <treatment_file>.contactloss.values.txt
#   (single column: 100 - col5 from processed/<new_name>/contacts_unisex.info)

set -euo pipefail

module purge
module load system
module load python/3.9.0
module load py-numpy/1.24.2_py39
module load py-pandas/2.0.1_py39

if [[ $# -ne 1 ]]; then
  echo "[ERROR] Usage: $0 <treatment_file.txt>" >&2
  exit 1
fi

TREAT_FILE="$1"
[[ -f "$TREAT_FILE" ]] || { echo "[ERROR] Not found: $TREAT_FILE" >&2; exit 1; }

BASE_DIR="$(cd "$(dirname "$TREAT_FILE")" && pwd)"
OUT_VALUES="${TREAT_FILE%.txt}.contactloss.values.txt"

# Find rename file in BASE_DIR (exactly one)
shopt -s nullglob
rename_files=( "${BASE_DIR}"/rename*.txt "${BASE_DIR}"/rename_*.txt )
shopt -u nullglob

if (( ${#rename_files[@]} > 1 )); then
  mapfile -t rename_files < <(printf "%s\n" "${rename_files[@]}" | awk '!seen[$0]++')
fi

if (( ${#rename_files[@]} != 1 )); then
  echo "[ERROR] Expected exactly 1 rename*.txt in: ${BASE_DIR}" >&2
  printf "  Found:\n" >&2
  printf "  %s\n" "${rename_files[@]:-NONE}" >&2
  exit 1
fi

RENAME_MAP="${rename_files[0]}"
PROCESSED_DIR="$(pwd)/processed"
INFO_FN="contacts_unisex.info"

echo "[INFO] TREAT_FILE   = ${TREAT_FILE}"
echo "[INFO] RENAME_MAP   = ${RENAME_MAP}"
echo "[INFO] PROCESSEDDIR = ${PROCESSED_DIR}"
echo "[INFO] OUTPUT_VALS  = ${OUT_VALUES}"

python3 - <<'PY' "$TREAT_FILE" "$RENAME_MAP" "$PROCESSED_DIR" "$INFO_FN" "$OUT_VALUES"
import sys, os

treat_file, rename_map, processed_dir, info_fn, out_vals = sys.argv[1:]

# load rename map
mp = {}
with open(rename_map, "r") as f:
    for ln in f:
        ln = ln.strip()
        if not ln:
            continue
        parts = ln.split()
        if len(parts) >= 2:
            mp[parts[0]] = parts[1]

# read treatment list and compute values
values = []
with open(treat_file, "r") as f:
    for ln in f:
        old = ln.strip()
        if not old:
            continue

        if old not in mp:
            raise KeyError(f"Not in rename map: {old}")

        nm = mp[old]
        info_path = os.path.join(processed_dir, nm, info_fn)

        if (not os.path.exists(info_path)) or os.path.getsize(info_path) == 0:
            raise FileNotFoundError(f"Missing/empty: {info_path}")

        with open(info_path, "r") as inf:
            line = inf.readline().strip()

        parts = line.split()
        if len(parts) < 5:
            raise ValueError(f"<5 columns in {info_path}: {line}")

        col5 = float(parts[4])
        values.append(100.0 - col5)

# write single-column file
os.makedirs(os.path.dirname(out_vals) or ".", exist_ok=True)
with open(out_vals, "w") as out:
    for v in values:
        out.write(f"{v:.6g}\n")

print(f"[DONE] wrote: {out_vals} ({len(values)} values)")
PY

echo "[DONE] Finished."

