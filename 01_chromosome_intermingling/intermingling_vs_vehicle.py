#!/usr/bin/env python3
"""
Chromosome intermingling: percentage of interchromosomal contacts per replicate,
compared between each treatment and vehicle (two-sided, equal-variance t-test,
BH correction per experiment).

Expected layout (one folder per experiment):
  <base>/experiment-XX_.../*_treatment_<name>.txt   one sample directory per line
  <sample_dir>/contacts_unisex.info                 hickit/dip-c contact summary;
                                                    field 5 = % intrachromosomal contacts
Exactly one treatment file per experiment must contain "vehicle" in its name.

Usage:
  python intermingling_vs_vehicle.py [base_dir] [output.tsv]
  (defaults: aux_data  contact_loss_summary.tsv)
"""
import pathlib
import sys
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

BASE = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "aux_data")

OUT_TSV = sys.argv[2] if len(sys.argv) > 2 else "contact_loss_summary.tsv"

results = []

def read_value(sample_dir: pathlib.Path):
    info = sample_dir / "contacts_unisex.info"
    if not info.exists():
        return None
    with open(info) as f:
        line = f.readline().strip().split()
    if len(line) < 5:
        return None
    try:
        v = float(line[4])
    except ValueError:
        return None
    return 100.0 - v


# ------------------------------------------------------------
# Iterate experiments
# ------------------------------------------------------------
for exp_dir in sorted(BASE.glob("experiment-*")):
    if not exp_dir.is_dir():
        continue

    exp_name = exp_dir.name

    # collect all treatment files
    treat_files = list(exp_dir.glob("*_treatment_*.txt"))
    if not treat_files:
        continue

    # identify vehicle
    vehicle_files = [f for f in treat_files if "vehicle" in f.name.lower()]
    if len(vehicle_files) != 1:
        print(f"[WARN] {exp_name}: expected 1 vehicle file, found {len(vehicle_files)}")
        continue

    vehicle_file = vehicle_files[0]

    # --------------------------------------------------------
    # Load vehicle values
    # --------------------------------------------------------
    vehicle_vals = []
    for line in vehicle_file.read_text().splitlines():
        p = pathlib.Path(line.strip())
        if not p.is_dir():
            continue
        v = read_value(p)
        if v is not None:
            vehicle_vals.append(v)

    vehicle_vals = np.array(vehicle_vals, dtype=float)

    if len(vehicle_vals) < 2:
        print(f"[WARN] {exp_name}: too few vehicle samples")
        continue

    # --------------------------------------------------------
    # Process each treatment
    # --------------------------------------------------------
    for tf in treat_files:
        if tf == vehicle_file:
            continue

        drug = tf.name.split("_treatment_", 1)[1].rsplit(".txt", 1)[0]

        vals = []
        for line in tf.read_text().splitlines():
            p = pathlib.Path(line.strip())
            if not p.is_dir():
                continue
            v = read_value(p)
            if v is not None:
                vals.append(v)

        vals = np.array(vals, dtype=float)

        if len(vals) < 2:
            continue

        # stats
        mean_val = vals.mean()
        n_treat = len(vals)
        n_vehicle = len(vehicle_vals)

        tstat, pval = ttest_ind(vals, vehicle_vals, equal_var=True)  # two-sided, equal-variance t-test

        results.append({
            "experiment": exp_name,
            "treatment": drug,
            "mean_100_minus_col5": mean_val,
            "n_treatment": n_treat,
            "n_vehicle": n_vehicle,
            "tstat": tstat,
            "pval": pval,
        })


# ------------------------------------------------------------
# Build dataframe + BH-FDR per experiment
# ------------------------------------------------------------
df = pd.DataFrame(results)

if df.empty:
    raise RuntimeError("No results generated.")

df["fdr"] = np.nan

for exp, g in df.groupby("experiment"):
    p = g["pval"].values
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(1, len(p) + 1))
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    df.loc[g.index[order], "fdr"] = np.clip(ranked, 0, 1)

df = df.sort_values(["experiment", "fdr", "pval"])

df.to_csv(OUT_TSV, sep="\t", index=False)
print(f"Wrote {OUT_TSV}")
