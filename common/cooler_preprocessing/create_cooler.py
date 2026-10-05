# import pandas as pd
# import cooler
# import sys
#
#!/usr/bin/env python3

import cooler
import pandas as pd
import numpy as np
import h5py
import sys

# ---------------- arguments ----------------
in_uri   = sys.argv[1]          # e.g. your.mcool::resolutions/5000
out_cool = sys.argv[2]
canonical = {str(c) for c in range(1, 23)} | {f"chr{c}" for c in range(1, 23)}
# -------------------------------------------

clr  = cooler.Cooler(in_uri)
bins = clr.bins()[:]

# 1. keep only wanted bins (canonical autosomes only)
mask        = bins["chrom"].isin(canonical)
old_idx     = bins[mask].index  # original bin IDs before reset
bins_keep   = bins[mask].reset_index(drop=True)
old2new     = pd.Series(np.arange(len(bins_keep)), index=old_idx)

print(f"kept {len(bins_keep)} / {len(bins)} bins")

# 2. load all pixels, filter, remap, and aggregate (non-chunked)
print("loading pixels...")
pix = clr.pixels()[:]
print(f"loaded {len(pix):,} pixels")

# filter: keep only pixels where both bins are in wanted chromosomes
m = mask.values[pix["bin1_id"]] & mask.values[pix["bin2_id"]]
pix = pix.loc[m].copy()
print(f"filtered to {len(pix):,} pixels")

# remap bin IDs to new indices
pix["bin1_id"] = pix["bin1_id"].map(old2new)
pix["bin2_id"] = pix["bin2_id"].map(old2new)

# sort by bin1_id, bin2_id (required by cooler)
pix.sort_values(["bin1_id", "bin2_id"], inplace=True)
pix.reset_index(drop=True, inplace=True)


# 3. write new cooler file
print(f"writing {out_cool}")

# get bin size (infer from bins since input cooler may not have it)
binsize = int(bins_keep["end"].iloc[0] - bins_keep["start"].iloc[0])

cooler.create_cooler(
    out_cool,
    bins=bins_keep,
    pixels=pix,
    ordered=True,
    ensure_sorted=True,
    dtypes={'count': np.int32},
)

# add bin-size attribute (required by cooler zoomify)
with h5py.File(out_cool, 'r+') as f:
    f.attrs['bin-size'] = binsize

print("done.")

# --- OLD CHUNKED APPROACH (commented out) ---
# block    = 1_000_000            # pixels per chunk
#
# # generator that yields remapped pixel chunks
# def pixel_chunks():
#     nnz = clr.info["nnz"]
#     for start in range(0, nnz, block):
#         stop  = min(start + block, nnz)
#         pix   = clr.pixels()[start:stop]
#         m     = mask.values[pix["bin1_id"]] & mask.values[pix["bin2_id"]]
#         if not m.any():
#             continue
#         pix   = pix.loc[m].copy()
#         pix["bin1_id"] = pix["bin1_id"].map(old2new)
#         pix["bin2_id"] = pix["bin2_id"].map(old2new)
#         # Aggregate duplicates within chunk by summing counts
#         pix = pix.groupby(["bin1_id", "bin2_id"], as_index=False)["count"].sum()
#         pix.sort_values(["bin1_id", "bin2_id"], inplace=True)
#         yield pix
#
# cooler.create_cooler(
#     out_cool,
#     bins=bins_keep,
#     pixels=pixel_chunks(),     # generator, low-RAM
#     ordered=True,
#     ensure_sorted=True,
#     dtypes={'count': np.int32},
# )








# bins_file = sys.argv[1]
# pixels_file = sys.argv[2]
# output_cool = sys.argv[3]
#
# # Read bins
# bins = pd.read_csv(bins_file, sep="\t", header=None, names=["chrom", "start", "end"])
# bins.reset_index(drop=True, inplace=True)
# n_bins = bins.shape[0]
#
# # Read pixels
# pixels = pd.read_csv(pixels_file, sep="\t", header=None, names=["bin1_id", "bin2_id", "count"])
#
# # Remap bin1_id and bin2_id to match 0-based new bin table
# # Assumes original bin1_id/bin2_id were from a larger table and now need shrinking
# # Create a mapping from old bin IDs (from full .cool) to new ones
# # Only needed if your pixel file is based on the old full bin index (very likely)
#
# # Build the mapping: old index → new index
# # This only works if you still know the original bin IDs before filtering
# # But since we don't — instead, just check if IDs are out of range
# if pixels[['bin1_id', 'bin2_id']].max().max() >= n_bins:
#     raise ValueError("Pixel file contains bin IDs beyond bin count. You must reindex the bins when filtering.")
#
# # If pixel bin IDs were already remapped, you're good
# cooler.create_cooler(
#     output_cool,
#     bins=bins,
#     pixels=pixels,
#     ordered=True
# )
