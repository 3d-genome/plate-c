#!/bin/bash
set -euo pipefail

module --force purge
module load devel
module load system
module load biology
module load math
module load python/3.9.0
module load py-numpy/1.24.2_py39
module load py-pandas/2.0.1_py39

if [[ $# -ne 1 ]]; then
  echo "[ERROR] Usage: $0 <cluster_file.txt>" >&2
  exit 1
fi

CLUSTER_FILE="$1"
[[ -f "$CLUSTER_FILE" ]] || { echo "[ERROR] Not found: $CLUSTER_FILE" >&2; exit 1; }

# -------------------------
# Config
# -------------------------
PROCESSED_BASE="/home/users/tttt/research/dip-c/processed"

# histogram parameters
min_distance=1000
binw=0.05
ld_min=2.95
ld_max=8.35

OUT_TSV="${CLUSTER_FILE%.txt}.log10_histogram_matrix.tsv"

echo "[INFO] CLUSTER_FILE   = ${CLUSTER_FILE}"
echo "[INFO] PROCESSED_BASE = ${PROCESSED_BASE}"
echo "[INFO] OUTPUT        = ${OUT_TSV}"

tmp_dir="$(mktemp -d)"
tmp_list="$(mktemp)"
tmp_missing="$(mktemp)"
trap 'rm -rf "$tmp_dir" "$tmp_list" "$tmp_missing"' EXIT

# -------------------------
# Step 1: build one temporary histogram per sample
# -------------------------
while IFS= read -r sample; do
  [[ -z "$sample" ]] && continue

  pairs_gz="${PROCESSED_BASE}/${sample}/contacts_unisex.pairs.gz"
  tmp_hist="${tmp_dir}/${sample}.hist.txt"

  if [[ ! -f "${pairs_gz}" ]]; then
    echo "[WARN] Missing input file: ${pairs_gz}" >> "$tmp_missing"
    continue
  fi

  echo "[INFO] Processing ${sample}"

  gunzip -c "$pairs_gz" \
    | grep -v "^#" \
    | awk -F $'\t' \
        -v min_distance="${min_distance}" \
        -v binw="${binw}" \
        -v ld_min="${ld_min}" \
        -v ld_max="${ld_max}" '
        BEGIN { total = 0 }

        ($2!="X"&&$2!="Y"&&$2!="chrX"&&$2!="chrY"&&
         $4!="X"&&$4!="Y"&&$4!="chrX"&&$4!="chrY") {

          # cis only
          if ($2 != $4) next

          d = $5 - $3
          if (d < 0) d = -d

          if (d >= min_distance) {
            ld = log(d)/log(10)
            if (ld < ld_min || ld > ld_max) next

            b = binw * int(ld / binw)
            hist[b]++
            total++
          }
        }

        END {
          if (total == 0) exit
          for (b in hist) {
            printf "%g %.8g\n", b, hist[b]/total
          }
        }
      ' \
    | sort -n > "$tmp_hist"

  if [[ -s "$tmp_hist" ]]; then
    printf "%s\t%s\n" "$sample" "$tmp_hist" >> "$tmp_list"
  else
    echo "[WARN] Empty histogram: ${sample}" >> "$tmp_missing"
    rm -f "$tmp_hist"
  fi

done < "$CLUSTER_FILE"

if [[ ! -s "$tmp_list" ]]; then
  echo "[ERROR] No valid histograms were generated." >&2
  [[ -s "$tmp_missing" ]] && cat "$tmp_missing" >&2
  exit 1
fi

# -------------------------
# Step 2: concatenate histograms into one matrix
# -------------------------
python3 - <<'PY' "$tmp_list" "$OUT_TSV"
import sys
import os
import pandas as pd
from functools import reduce

lst_path = sys.argv[1]
out_tsv  = sys.argv[2]

items = []
with open(lst_path, "r") as f:
    for ln in f:
        ln = ln.strip()
        if not ln:
            continue
        sample, path = ln.split("\t", 1)
        items.append((sample, path))

dfs = []
for sample, path in items:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        continue
    df = pd.read_csv(path, sep=r"\s+", header=None, names=["x", sample])
    df["x"] = pd.to_numeric(df["x"], errors="coerce")
    df[sample] = pd.to_numeric(df[sample], errors="coerce")
    df = df.dropna(subset=["x"])
    dfs.append(df)

if not dfs:
    sys.stderr.write("[ERROR] No readable temporary histograms.\n")
    sys.exit(2)

merged = reduce(lambda l, r: pd.merge(l, r, on="x", how="outer"), dfs)
merged = merged.sort_values("x")
merged.to_csv(out_tsv, sep="\t", index=False)
PY

if [[ -s "$tmp_missing" ]]; then
  echo "[WARN] Some samples were skipped:" >&2
  cat "$tmp_missing" >&2
fi

echo "[DONE] Wrote matrix: $OUT_TSV"
