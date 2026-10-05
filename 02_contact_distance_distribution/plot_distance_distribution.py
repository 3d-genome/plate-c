#!/usr/bin/env python3
"""
Plot contact-distance distributions (log10 genomic distance) for vehicle and treatments.

Input: one "*.log10_histogram_matrix.tsv" per treatment. Column 1 holds the
log10(distance in bp) bin centres ("x"); every further column is one
technical replicate (normalized contact frequency per bin).

Each replicate is drawn as a faint line and the replicate mean as a solid line,
as in the contact-distance panels of Figs. 3e, 5b and 7f.

Example:
  python 02_contact_distance_distribution/plot_distance_distribution.py \
      demo/data/distance_histograms/*DMSO*.tsv demo/data/distance_histograms/*SGI-1027*.tsv \
      --colors DMSO=black SGI-1027=#dc0000 --output distance_distribution.svg
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def drug_from_filename(path):
    bn = os.path.basename(path)
    return bn.split("_treatment_")[1].split(".log10")[0] if "_treatment_" in bn else bn


def load_matrix(path, x_min, x_max):
    df = pd.read_csv(path, sep="\t", header=None)
    x = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy()
    Y = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy()
    ok = np.isfinite(x) & (x >= 0)  # drops the header row
    x, Y = x[ok], Y[ok, :]
    sel = (x >= x_min) & (x <= x_max)
    return x[sel], Y[sel, :]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--colors", nargs="*", default=[], help="DRUG=COLOR pairs")
    ap.add_argument("--x-min", type=float, default=2.95)
    ap.add_argument("--x-max", type=float, default=8.35)
    ap.add_argument("--output", default="distance_distribution.svg")
    ap.add_argument("--summary", help="Optional TSV with the mean curve per treatment")
    a = ap.parse_args(argv)

    colors = dict(c.split("=", 1) for c in a.colors)
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    summary = {}
    for f in a.files:
        drug = drug_from_filename(f)
        color = colors.get(drug, "black")
        x, Y = load_matrix(f, a.x_min, a.x_max)
        for j in range(Y.shape[1]):
            ax.plot(x, Y[:, j], color=color, alpha=0.15, linewidth=0.6, zorder=1)
        mean = np.nanmean(Y, axis=1)
        ax.plot(x, mean, color=color, linewidth=1.8, zorder=2, label=f"{drug} (n={Y.shape[1]})")
        summary["log10_distance"] = x
        summary[f"{drug}_mean"] = mean

    ax.set_xlim(a.x_min, a.x_max)
    ax.set_xlabel("log10 contact distance (bp)")
    ax.set_ylabel("Contact frequency")
    ax.legend(frameon=False, fontsize=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(a.output, bbox_inches="tight")
    print(f"Wrote {a.output}")
    if a.summary:
        pd.DataFrame(summary).to_csv(a.summary, sep="\t", index=False, float_format="%.6g")
        print(f"Wrote {a.summary}")


if __name__ == "__main__":
    main()
