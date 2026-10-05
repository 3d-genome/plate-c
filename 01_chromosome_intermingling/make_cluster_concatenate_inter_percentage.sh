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

PROCESSED_DIR="/home/users/tttt/research/dip-c/processed"
INFO_FN="contacts_unisex.info"

OUT_VALUES="${CLUSTER_FILE%.txt}.inter_percentage.txt"

echo "[INFO] CLUSTER_FILE = ${CLUSTER_FILE}"
echo "[INFO] PROCESSEDDIR = ${PROCESSED_DIR}"
echo "[INFO] OUTPUT_VALS  = ${OUT_VALUES}"

python3 - <<'PY' "$CLUSTER_FILE" "$PROCESSED_DIR" "$INFO_FN" "$OUT_VALUES"
import sys, os

cluster_file, processed_dir, info_fn, out_vals = sys.argv[1:]

values = []

with open(cluster_file) as f:
    for ln in f:
        sample = ln.strip()
        if not sample:
            continue

        info_path = os.path.join(processed_dir, sample, info_fn)

        if (not os.path.exists(info_path)) or os.path.getsize(info_path) == 0:
            raise FileNotFoundError(f"Missing/empty: {info_path}")

        with open(info_path) as inf:
            line = inf.readline().strip()

        parts = line.split()

        if len(parts) < 5:
            raise ValueError(f"<5 columns in {info_path}: {line}")

        col5 = float(parts[4])
        values.append(100.0 - col5)

os.makedirs(os.path.dirname(out_vals) or ".", exist_ok=True)

with open(out_vals, "w") as out:
    for v in values:
        out.write(f"{v:.6g}\n")

print(f"[DONE] wrote: {out_vals} ({len(values)} values)")
PY
