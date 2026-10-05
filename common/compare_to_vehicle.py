#!/usr/bin/env python3
"""
Compare a per-replicate architectural attribute between each treatment and vehicle.

This is the shared statistics step used for the attribute folders 01-04, 06 and 07:
for every treatment, its technical replicates are compared with the vehicle
replicates of the same experiment by a two-sided t-test, and p-values are
corrected across all treatments of the experiment with the Benjamini-Hochberg
procedure.

Input: a tab-separated table with one row per replicate (sample), containing
  - a treatment column (e.g. "drug" or "treatment_index"),
  - one or more numeric attribute columns (e.g. "saddle_strength_extent10_log2").

Output: one row per (treatment, attribute) with n, mean, s.e.m., p-value,
BH-adjusted q-value, and -log10(q).

Example:
  python common/compare_to_vehicle.py \
      --input demo/data/saddle_strength_example.tsv \
      --treatment-col drug --vehicle DMSO \
      --value-cols saddle_strength_extent10_log2 relative_strength_extent10_log2 \
      --output saddle_stats.tsv

The test is a two-sided, equal-variance (Student's) t-test, as stated in the
figure legends. --welch switches to Welch's unequal-variance t-test.
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values; NaNs are passed through."""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = q.size
    if n == 0:
        return out
    order = np.argsort(q)
    ranked = q[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def compare(df, treatment_col, vehicle_labels, value_cols, equal_var=True, min_n=2):
    is_vehicle = df[treatment_col].astype(str).isin(vehicle_labels)
    vehicle = df[is_vehicle]
    if len(vehicle) < min_n:
        raise ValueError(f"Found {len(vehicle)} vehicle replicates for labels {vehicle_labels}; need >= {min_n}.")

    rows = []
    for col in value_cols:
        v = pd.to_numeric(vehicle[col], errors="coerce").dropna().to_numpy()
        for treatment, g in df[~is_vehicle].groupby(treatment_col, sort=True):
            x = pd.to_numeric(g[col], errors="coerce").dropna().to_numpy()
            p = ttest_ind(x, v, equal_var=equal_var).pvalue if x.size >= min_n else np.nan
            rows.append({
                "attribute": col,
                "treatment": treatment,
                "n": x.size,
                "mean": x.mean() if x.size else np.nan,
                "sem": x.std(ddof=1) / np.sqrt(x.size) if x.size > 1 else np.nan,
                "n_vehicle": v.size,
                "mean_vehicle": v.mean(),
                "delta_vs_vehicle": (x.mean() - v.mean()) if x.size else np.nan,
                "p_value": p,
            })

    res = pd.DataFrame(rows)
    # BH correction across treatments, separately for each attribute
    res["q_value"] = np.nan
    for col, idx in res.groupby("attribute").groups.items():
        res.loc[idx, "q_value"] = bh_fdr(res.loc[idx, "p_value"].to_numpy())
    res["neg_log10_q"] = -np.log10(res["q_value"])
    return res.sort_values(["attribute", "p_value"]).reset_index(drop=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="Per-replicate TSV")
    ap.add_argument("--treatment-col", required=True)
    ap.add_argument("--vehicle", nargs="+", default=["DMSO"], help="Label(s) of vehicle replicates")
    ap.add_argument("--value-cols", nargs="+", required=True)
    ap.add_argument("--welch", action="store_true", help="Use Welch's unequal-variance t-test instead of the default equal-variance t-test")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)

    df = pd.read_csv(a.input, sep="\t")
    missing = [c for c in [a.treatment_col] + a.value_cols if c not in df.columns]
    if missing:
        sys.exit(f"Missing columns in {a.input}: {missing}")
    res = compare(df, a.treatment_col, a.vehicle, a.value_cols, equal_var=not a.welch)
    res.to_csv(a.output, sep="\t", index=False, float_format="%.6g")
    print(f"Wrote {a.output} ({len(res)} rows)")


if __name__ == "__main__":
    main()
