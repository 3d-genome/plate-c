#!/bin/bash
# concat_treatment_log10_histogram.sh
#
# Usage:
#   ./concat_treatment_log10_histogram.sh aux_data/<experiment_dir>/<treatment_file>.txt
#
# What it does:
# - Reads a treatment file that contains old ABS paths (one per line)
# - Finds rename_*.txt in the same directory as the treatment file (BASE)
# - Uses rename map: old_abs_path -> new_name
# - For each replicate, loads:
#     processed/<new_name>/contacts_unisex.distance_log10_cis_histogram.txt
#   which is assumed to be 2 columns: (bin_x, value)
# - Outputs a matrix TSV with:
#     bin_x    <new_name_1>  <new_name_2>  ...
#   aligned by bin_x (outer join across files)

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
TREAT_BN="$(basename "$TREAT_FILE")"
OUT_TSV="${TREAT_FILE%.txt}.log10_histogram_matrix.tsv"

# Find rename file in BASE_DIR (exactly one)
shopt -s nullglob
rename_files=( "${BASE_DIR}"/rename*.txt "${BASE_DIR}"/rename_*.txt )
shopt -u nullglob

# de-dup (in case both globs hit the same file)
if (( ${#rename_files[@]} > 1 )); then
  # keep unique
  mapfile -t rename_files < <(printf "%s\n" "${rename_files[@]}" | awk '!seen[$0]++')
fi

if (( ${#rename_files[@]} != 1 )); then
  echo "[ERROR] Expected exactly 1 rename*.txt in: ${BASE_DIR}" >&2
  printf "  Found:\n" >&2
  printf "  %s\n" "${rename_files[@]:-NONE}" >&2
  exit 1
fi

RENAME_MAP="${rename_files[0]}"

# You said histograms live here:
PROCESSED_DIR="$(pwd)/processed"
HIST_FN="contacts_unisex.distance_log10_cis_histogram.txt"

echo "[INFO] TREAT_FILE   = ${TREAT_FILE}"
echo "[INFO] BASE_DIR     = ${BASE_DIR}"
echo "[INFO] RENAME_MAP   = ${RENAME_MAP}"
echo "[INFO] PROCESSEDDIR = ${PROCESSED_DIR}"
echo "[INFO] OUTPUT       = ${OUT_TSV}"

# Build a list of histogram paths + column names in a temp file
tmp_list="$(mktemp)"
tmp_warn="$(mktemp)"
trap 'rm -f "$tmp_list" "$tmp_warn"' EXIT

# Parse rename map into a fast lookup in awk, then map each old path -> new_name
# and emit: new_name<TAB>hist_path
awk -v RENAME_MAP="$RENAME_MAP" \
    -v PROCESSED_DIR="$PROCESSED_DIR" \
    -v HIST_FN="$HIST_FN" \
    -v WARN_FILE="$tmp_warn" '
BEGIN{
  # load rename map: old -> new
  while ((getline line < RENAME_MAP) > 0) {
    gsub(/\r$/, "", line)
    if (line == "") continue
    n = split(line, a, /[ \t]+/)
    if (n >= 2) mp[a[1]] = a[2]
  }
  close(RENAME_MAP)
}
{
  gsub(/\r$/, "", $0)
  if ($0 == "") next
  old = $0
  if (!(old in mp)) {
    print "[WARN] Not in rename map: " old > WARN_FILE
    next
  }
  nm = mp[old]
  hist = PROCESSED_DIR "/" nm "/" HIST_FN
  # We do not check existence here (python will), but we can still warn early:
  print nm "\t" hist
}
' "$TREAT_FILE" > "$tmp_list"

if [[ ! -s "$tmp_list" ]]; then
  echo "[ERROR] No valid entries after mapping old paths -> new_name (check rename map + treatment file)." >&2
  [[ -s "$tmp_warn" ]] && cat "$tmp_warn" >&2
  exit 1
fi

# Do the alignment/join in python (robust outer-join on first column)
python3 - <<'PY' "$tmp_list" "$OUT_TSV"
import sys, os
import pandas as pd

lst_path = sys.argv[1]
out_tsv  = sys.argv[2]

items = []
with open(lst_path, "r") as f:
    for ln in f:
        ln = ln.strip()
        if not ln:
            continue
        nm, path = ln.split("\t", 1)
        items.append((nm, path))

dfs = []
missing = []
for nm, path in items:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        missing.append((nm, path))
        continue
    # histogram file: two columns (x, y). allow whitespace delimiter.
    df = pd.read_csv(path, sep=r"\s+", header=None, names=["x", nm])
    # ensure numeric, drop bad lines
    df["x"] = pd.to_numeric(df["x"], errors="coerce")
    df[nm]  = pd.to_numeric(df[nm], errors="coerce")
    df = df.dropna(subset=["x"])
    dfs.append(df)

if not dfs:
    sys.stderr.write("[ERROR] No histogram files found/readable.\n")
    for nm, path in missing:
        sys.stderr.write(f"  missing/empty: {nm}\t{path}\n")
    sys.exit(2)

# Outer-merge on x
merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df, on="x", how="outer")

merged = merged.sort_values("x")

# Write
merged.to_csv(out_tsv, sep="\t", index=False)

# Emit warnings
if missing:
    sys.stderr.write("[WARN] Some histograms were missing/empty and were skipped:\n")
    for nm, path in missing:
        sys.stderr.write(f"  {nm}\t{path}\n")
PY

# Print any rename-map warnings
if [[ -s "$tmp_warn" ]]; then
  echo "[WARN] Issues while mapping old paths -> new_name:" >&2
  cat "$tmp_warn" >&2
fi

echo "[DONE] Wrote matrix: $OUT_TSV"

