import warnings
warnings.filterwarnings("ignore")
from itertools import combinations
import os
import sys
from pathlib import Path
import re
import logging
logging.basicConfig(level=logging.INFO)
from functools import partial

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import gridspec
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
from matplotlib.ticker import EngFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.filters import threshold_li, threshold_otsu # type: ignore
from matplotlib.patches import Rectangle

import numpy as np
import pandas as pd
from multiprocessing import Pool
import bioframe # type: ignore
import cooler # type: ignore
import cooltools # type: ignore
from cooltools.lib.numutils import fill_diag  # type: ignore
from cooltools.lib._query import CSRSelector # type: ignore
from cooltools.lib import peaks, numutils # type: ignore
from cooltools.lib.checks import is_compatible_viewframe, is_cooler_balanced # type: ignore
from cooltools.lib.common import make_cooler_view, pool_decorator # type: ignore

from coolpuppy import coolpup # type: ignore
from coolpuppy import plotpup # pyright: ignore[reportMissingImports]
import random
from packaging import version
bp_formatter = EngFormatter('b')

if version.parse(cooltools.__version__) < version.parse('0.5.2'):
    raise AssertionError("tutorial relies on cooltools version 0.5.2 or higher,"+
                         "please check your cooltools version and update to the latest")

def aggregate_plot1(data, cmap, resolution, output_dir, title1, title2, vmax = None):

        stack = np.stack(data)     # (N, H, W)
        mean_matrix = np.nanmean(stack, axis=0)
        log_map = np.log2(mean_matrix)
        # Plot setup
        fig = plt.figure(figsize=(5, 6))
        gs = gridspec.GridSpec(2, 1, height_ratios=[20, 1], hspace=.6)

        ax = fig.add_subplot(gs[0])
        cbar_ax = fig.add_subplot(gs[1])

        # Show image on main ax
        im = ax.imshow(log_map,
                       cmap=cmap,
                       interpolation='none',
                       vmin=-.75, vmax=.75
                       )

        # Colorbar underneath
        cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar.set_label(r"$log_{2}(mean O/E)$")
        # Distance-based tick labels
        # center = mean_matrix.shape[0] // 2
        # ticks_pixels = np.linspace(0, mean_matrix.shape[0]-1, 5)
        # ticks_kb = ((ticks_pixels - center) * resolution)//1000

        n   = mean_matrix.shape[0]          # e.g. 200
        mid = n // 2
        ticks_px   = np.linspace(0, n-1, 5)       # pixel positions
        ticks_frac = (ticks_px - mid) / mid        # –1 … +1
        ticks_pct  = ticks_frac * 100              # –100 … +100 %

        ax.set_xticks(ticks_px)
        ax.set_yticks(ticks_px)
        ax.set_xticklabels(np.round(ticks_frac, 2))
        ax.set_yticklabels(np.round(ticks_frac, 2))
        ax.set_xlabel("Relative position (boundary = ±1)", fontsize = 12)
        ax.set_ylabel("Relative position (boundary = ±1)", fontsize = 12)

        # ax.set_xticks(ticks_pixels)
        # ax.set_xticklabels(ticks_kb.astype(int))
        # ax.set_yticks(ticks_pixels)
        # ax.set_yticklabels(ticks_kb.astype(int))
        # ax.set_xlabel("Distance from boundary (kb)")
        # ax.set_ylabel("Distance from boundary (kb)")

        ax.set_title(title1 + " " + str(resolution) + "kb Resolution: \n" +  title2 + " Identification and Aggregation", fontsize=12)
        plt.tight_layout()
        save_path = output_dir+'/TAD_' + title1 + '_' + str(resolution) + '_aggrplot' + '.png'
        save_path = save_path.replace('TADs', title2)
        plt.savefig(save_path)

def aggregate_plot2(data, cmap, resolution, output_dir, title1, title2):

        # vmax_ref = 3.5       # colour limit you liked at 50 kb
        # ref_res  = 25_000      # reference resolution in bp
        #
        # scale = (resolution / ref_res)**2
        # vmax  = vmax_ref * scale
        # vmin  = -vmax          # keep symmetric colour‑bar

        stack = np.stack(data)     # (N, H, W)
        mean_matrix = np.nanmean(stack, axis=0)
        log_map = np.log2(mean_matrix)
        # Plot setup
        fig = plt.figure(figsize=(5, 6))
        gs = gridspec.GridSpec(2, 1, height_ratios=[20, 1], hspace=.6)

        ax = fig.add_subplot(gs[0])
        cbar_ax = fig.add_subplot(gs[1])

        # Show image on main ax
        im = ax.imshow(log_map,
                       cmap=cmap,
                       interpolation='none',
                       vmin=0, vmax=5
                       )

        # Colorbar underneath
        cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar.set_label(r"$log_{2}(mean O/E)$")
        # Distance-based tick labels
        # center = mean_matrix.shape[0] // 2
        # ticks_pixels = np.linspace(0, mean_matrix.shape[0]-1, 5)
        # ticks_kb = ((ticks_pixels - center) * resolution)//1000

        n   = mean_matrix.shape[0]          # e.g. 200
        mid = n // 2
        ticks_px   = np.linspace(0, n-1, 5)       # pixel positions
        ticks_frac = (ticks_px - mid) / mid        # –1 … +1
        ticks_pct  = ticks_frac * 100              # –100 … +100 %

        ax.set_xticks(ticks_px)
        ax.set_yticks(ticks_px)
        ax.set_xticklabels(np.round(ticks_frac, 2))
        ax.set_yticklabels(np.round(ticks_frac, 2))
        ax.set_xlabel("Relative position (boundary = ±1)", fontsize = 12)
        ax.set_ylabel("Relative position (boundary = ±1)", fontsize = 12)

        # ax.set_xticks(ticks_pixels)
        # ax.set_xticklabels(ticks_kb.astype(int))
        # ax.set_yticks(ticks_pixels)
        # ax.set_yticklabels(ticks_kb.astype(int))
        # ax.set_xlabel("Distance from boundary (kb)")
        # ax.set_ylabel("Distance from boundary (kb)")

        ax.set_title(title1 + " " + str(resolution) + "kb Resolution: \n" +  title2 + " Identification and Aggregation", fontsize=12)
        plt.tight_layout()
        save_path = output_dir+'/Loop_' + title1 + '_' + str(resolution) + '_aggrplot' + '.png'
        save_path = save_path.replace('TADs', title2)
        plt.savefig(save_path)

def TAD_sanity_check(data, domains, chr, start, end, resolution, path):
  resolution = resolution
  region = (chr, start, end)

  df_region = domains[
      (domains['chrom'] == region[0]) &
      (domains['start'] >= region[1]) &
      (domains['end'] <= region[2])
  ].copy()

  norm = LogNorm(vmax=0.1, vmin=0.001)
  data = data.matrix(balance=True).fetch(region)

  f, ax = plt.subplots(figsize=(18, 10))
  im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
  plt.title(f'{chr}:{start}-{end}')
  ax.set_aspect(0.5)
  ax.set_ylim(0, 50*resolution)
  format_ticks(ax, rotate=False)
  ax.xaxis.set_visible(False)

  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
  plt.colorbar(im, cax=cax)

  idx = resolution
  offset  = region[1]
  max_pos = (region[2] - region[1]) // resolution      # 1500 bins for 15 Mb
  contact_matrix = np.full((max_pos, max_pos), np.nan) # start with NaNs

  for _, row in df_region[:idx].iterrows():
      start = (row['start'] - offset) // resolution
      end = (row['end'] - offset) // resolution

      if end - start > 2:
          contact_matrix[start:end, start:end] = 1
          contact_matrix[start+1:end-1, start+1:end-1] = np.nan
      else:
          contact_matrix[start:end, start:end] = 1

  im2 = pcolormesh_45deg(ax, contact_matrix, start=region[1], resolution=resolution, cmap='Blues', vmax=1, vmin=-1, alpha=0.6)
  plt.tight_layout()
  plt.savefig(path + '/' + f'{resolution}_{chr}_{start}_{end}' + '_random_TAD_sanitycheck_plot.pdf')

def loops_sanity_check(data, loops, chr, start, end, resolution, path):
  # define a region to look into as an example
  resolution = resolution
  region = (chr, start, end)

  # heatmap kwargs
  matshow_kwargs = dict(
      cmap='gist_yarg',
      norm=LogNorm(vmax=0.05),
      extent=(start, end, end, start)
  )

  # colorbar kwargs
  colorbar_kwargs = dict(fraction=0.046, label='corrected frequencies')

  # compute heatmap for the region
  region_matrix = data.matrix(balance=True).fetch(region)
  for diag in [-1,0,1]:
      region_matrix = fill_diag(region_matrix, np.nan, i=diag)

  # see viz.ipynb for details of heatmap visualization
  f, ax = plt.subplots(figsize=(7,7))
  im = ax.matshow( region_matrix, **matshow_kwargs)
  format_ticks(ax, rotate=False)
  plt.title(f'{chr}:{start}-{end}')
  plt.colorbar(im, ax=ax, **colorbar_kwargs)

  # draw rectangular "boxes" around pixels called as dots in the "region":
  for box in rectangles_around_dots(loops, region, lw=1.5):
      ax.add_patch(box)
  plt.savefig(path.replace('TAD', 'Loop') + '/' + f'{resolution}_{chr}_{start}_{end}' + '_random_loop_sanitycheck_plot.pdf')

def random_pair(high=1757):
    a = random.randint(0, high - 8)
    return a, a + 8

def rectangles_around_dots(dots_df, region, loc="upper", lw=1, ec="cyan", fc="none"):
    """
    yield a series of rectangles around called dots in a given region
    """
    # select dots from the region:
    df_reg = bioframe.select(
        bioframe.select(dots_df, region, cols=("chrom1","start1","end1")),
        region,
        cols=("chrom2","start2","end2"),
    )
    rectangle_kwargs = dict(lw=lw, ec=ec, fc=fc)
    # draw rectangular "boxes" around pixels called as dots in the "region":
    for s1, s2, e1, e2 in df_reg[["start1", "start2", "end1", "end2"]].itertuples(index=False):
        width1 = e1 - s1
        width2 = e2 - s2
        if loc == "upper":
            yield patches.Rectangle((s2, s1), width2, width1, **rectangle_kwargs)
        elif loc == "lower":
            yield patches.Rectangle((s1, s2), width1, width2, **rectangle_kwargs)
        else:
            raise ValueError("loc has to be uppper or lower")

def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    import itertools
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im

def format_ticks(ax, x=True, y=True, rotate=True):
    """format ticks with genomic coordinates as human readable"""
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

# def get_boundary_strength(clr, merged_bound: pd.DataFrame, flank: int = 300_000, ignore_diags: int = 2):
#
#     """
#     Snap boundary centres to the bin grid, recompute insulation at `flank`,
#     and return boundary-strength values for those reference boundaries.
#
#     Parameters
#     ----------
#     clr : cooler.Cooler
#         Opened Cooler object (any resolution).
#     merged_bound : pd.DataFrame
#         BED-like table with columns 'chrom', 'start' (and optionally 'end').
#         Coordinates may be bin centres.
#     flank : int, default 300_000
#         Window size (bp) used for insulation / boundary strength.
#     ignore_diags : int, default 2
#         Passed through to `cooltools.insulation`.
#
#     Returns
#     -------
#     pd.DataFrame
#         Columns: chrom, bin_start, boundary_strength_<flank>
#         Only rows that overlap the reference boundaries.
#     """
#     bin_size = clr.binsize                              # e.g. 25_000
#
#     # snap centres → nearest bin edge and deduplicate
#     bnds = (merged_bound
#             .assign(bin_start=((merged_bound.start + bin_size // 2) // bin_size) * bin_size)
#             .loc[:, ["chrom", "bin_start"]]
#             .drop_duplicates())
#
#     # recompute insulation / boundary strength
#     ins = (cooltools.insulation(
#                 clr,
#                 window_bp=[flank],
#                 ignore_diags=ignore_diags,
#                 append_raw_scores=True)
#             .rename(columns={"start": "bin_start"})
#             .loc[:, ["chrom", "bin_start", f"boundary_strength_{flank}"]])
#
#     # intersect with reference boundaries
#     merged = ins.merge(bnds, on=["chrom", "bin_start"], how="inner")
#     print(f"matches: {len(merged)} / {len(bnds)}")  # quick QC
#     return merged

# def get_boundary_strength(clr, merged_bound: pd.DataFrame, flank: int = 300_000, ignore_diags: int = 2):
#
#     """
#     Snap boundary centres to the bin grid, recompute insulation at `flank`,
#     and return boundary-strength values for those reference boundaries.
#
#     Parameters
#     ----------
#     clr : cooler.Cooler
#         Opened Cooler object (any resolution).
#     merged_bound : pd.DataFrame
#         BED-like table with columns 'chrom', 'start' (and optionally 'end').
#         Coordinates may be bin centres.
#     flank : int, default 300_000
#         Window size (bp) used for insulation / boundary strength.
#     ignore_diags : int, default 2
#         Passed through to `cooltools.insulation`.
#
#     Returns
#     -------
#     pd.DataFrame
#         Columns: chrom, bin_start, boundary_strength_<flank>
#         Only rows that overlap the reference boundaries.
#     """
#     bin_size = clr.binsize                              # e.g. 25_000
#
#     # snap centres → nearest bin edge and deduplicate
#     bnds = (merged_bound[merged_bound.score < -0.2]
#             .assign(bin_start=((merged_bound.start + bin_size // 2) // bin_size) * bin_size)
#             .loc[:, ["chrom", "bin_start"]]
#             .drop_duplicates())
#
#     w_bins = flank // bin_size       # number of bins in the diamond
#     diamond_area = w_bins * w_bins
#     # recompute insulation / boundary strength
#     ins = (cooltools.insulation(
#                 clr,
#                 window_bp=[flank],
#                 ignore_diags=ignore_diags,
#                 append_raw_scores=True))
#
#     ins = ins[ins.n_valid_pixels_300000 > 0.5 * diamond_area]
#     ins = ins.rename(columns={"start": "bin_start"}).loc[:, ["chrom", "bin_start", f"boundary_strength_{flank}"]]
#
#     # intersect with reference boundaries
#     merged = ins.merge(bnds, on=["chrom", "bin_start"], how="inner")
#     print(f"matches: {len(merged)} / {len(bnds)}")  # quick QC
#     return merged

# def boundary_strength_summary(clr, merged_bound, flank=300_000):
#     """Return one row of summary stats for a given cooler & BED file."""
#     bs      = get_boundary_strength(clr, merged_bound, flank)   # your old fn
#     col     = f"boundary_strength_{flank}"
#     # called  = bs.query("is_boundary")
#     return dict(
#         n_bound      = len(bs),
#         # n_called     = len(called),
#         # frac_ret     = len(called)/len(bs),
#         mean_bs_all  = bs[col].mean(),
#         median_bs_all= bs[col].median(),
#         # mean_bs_call = called[col].mean(),
#         # median_bs_call=called[col].median()
#     )


def get_boundary_strength(clr, merged_bound: pd.DataFrame, flank: int = 300_000, ignore_diags: int = 2):

    """
    Snap boundary centres to the bin grid, recompute insulation at `flank`,
    and return boundary-strength values for those reference boundaries.

    Parameters
    ----------
    clr : cooler.Cooler
        Opened Cooler object (any resolution).
    merged_bound : pd.DataFrame
        BED-like table with columns 'chrom', 'start' (and optionally 'end').
        Coordinates may be bin centres.
    flank : int, default 300_000
        Window size (bp) used for insulation / boundary strength.
    ignore_diags : int, default 2
        Passed through to `cooltools.insulation`.

    Returns
    -------
    pd.DataFrame
        Columns: chrom, bin_start, boundary_strength_<flank>
        Only rows that overlap the reference boundaries.
    """
    bin_size = clr.binsize                              # e.g. 25_000

    # snap centres → nearest bin edge and deduplicate
    bnds = (merged_bound #[merged_bound.score < -0.2]
            .assign(bin_start=((merged_bound.start + bin_size // 2) // bin_size) * bin_size)
            .loc[:, ["chrom", "bin_start"]]
            .drop_duplicates())

    windows = [flank]
    print(windows)
    w_bins = windows[-1] // bin_size
    diamond_area = w_bins * w_bins
    print(diamond_area)

    # recompute insulation / boundary strength
    ins0 = (insulation( #my version of insulation!
                clr,
                window_bp=[windows[-1]],
                ignore_diags=ignore_diags,
                append_raw_scores=True))

    ins = ins0[(ins0.is_bad_bin == False) & (ins0[f'n_valid_pixels_{windows[-1]}'] > 0.50 * diamond_area)]
    ins = ins.rename(columns={"start": "bin_start"}).loc[:, ["chrom", "bin_start", f'insulation_score_{windows[-1]}', f'n_valid_pixels_{windows[-1]}']]
    # intersect with reference boundaries
    merged = ins.merge(bnds, on=["chrom", "bin_start"], how="inner")
    return merged, ins0


def boundary_strength_summary(clr, merged_bound, path, flank=300000):
    """Return one row of summary stats for a given cooler & BED file."""
    merged, ins      = get_boundary_strength(clr, merged_bound, flank)
    col = f"insulation_score_{flank}"
    merged = merged.replace([np.inf, -np.inf], np.nan).dropna() #these are windows that were sparse
    print(merged[col].mean())
    print(merged[col].median())

    ins.to_csv(path + '_insulation0_table.csv')
    merged.to_csv(path + '_merged_insulation_table.csv')

    print(np.log2(merged[f"insulation_score_{flank}"].mean()))
    return dict(
        n_bound      = len(merged),
        # n_called     = len(called),
        # frac_ret     = len(called)/len(bs),
        log2_mean_aggr_ins_all  = np.log2(merged[col].mean()),
        log2_median_aggr_ins_all= np.log2(merged[col].median()))
        # mean_bs_call = called[col].mean(),
        # median_bs_call=called[col].median())

# def loop_strength_mean_map(mean_map):
#     print(mean_map)
#     mid   = mean_map.shape[0] // 2
#     peak  = mean_map[mid-1:mid+2, mid-1:mid+2].mean()
#     print('peak: ' + str(peak))
#     ll    = mean_map[mid+2:mid+4, mid-4:mid-2].mean()
#     print('ll: ' + str(ll))
#     return dict(P2LL_global = float(peak-ll))

# def loop_strength_summary(loops_matrix):
#     mats   = np.stack(loops_matrix)    # (N, H, W)
#     n      = mats.shape[1]; mid = n//2
#
#     peak = mats[:, mid-1:mid+2, mid-1:mid+2].mean(axis=(1,2))
#     ll   = mats[:, mid+2:mid+4,   mid-4:mid-2].mean(axis=(1,2))
#     p2ll = np.log2(peak/ll)
#
#     summary = dict(n_loops=len(mats),
#                   mean_P2LL=p2ll.mean(),
#                   median_P2LL=np.median(p2ll))
#     return summary


# --- 1) P2LL with pixel-fraction scaling ---
def p2ll_px(pup_df, peak_frac=1/4, bg_frac=1/4, off_frac=1/4,
            force_odd_peak=True, force_odd_bg=False, push_ll_to_corner=True):

    # stack = np.stack(mean_map['data'])     # (N, H, W)
    # mean_matrix = np.nanmean(mean_map['data'], axis=0)
    mean_map = np.log2(pup_df['data'][0])

    n = mean_map.shape[0]
    mid = (n - 1) // 2
    odd = (lambda w: w if w % 2 else w + 1)

    wP  = int(round(n * peak_frac)); wP = odd(max(1, wP)) if force_odd_peak else max(1, wP)
    wB  = int(round(n * bg_frac));   wB = odd(max(1, wB)) if force_odd_bg   else max(1, wB)
    off = max(1, int(round(n * off_frac)))

    # keep BG outside peak and inside the map
    need = (wP // 2) + (wB // 2) + 1
    off = max(off, need)
    off_max = min(n - 1 - (wB // 2) - mid,   # row: mid+off <= n-1-halfB
                  mid - (wB // 2))           # col: mid-off >= halfB
    off = max(1, min(off, off_max))

    # optionally snap LL to the true lower-left corner
    if push_ll_to_corner:
        hB = wB // 2
        off_corner = mid - hB
        off_max = min(n - 1 - hB - mid,  mid - hB)
        off = max(1, min(off_corner, off_max))

    h = lambda w: w // 2
    def box(rc, cc, w):
        r0, r1 = rc - h(w), rc + h(w) + 1
        c0, c1 = cc - h(w), cc + h(w) + 1
        return slice(max(0, r0), min(n, r1)), slice(max(0, c0), min(n, c1))

    rp, cp   = box(mid, mid, wP)                 # peak (center)
    rll, cll = box(mid + off, mid - off, wB)     # lower-left BG

    peak = mean_map[rp, cp].mean()
    ll   = mean_map[rll, cll].mean()
    out = {
        "P2LL": float(peak - ll),
        "peak": float(peak), "LL": float(ll),
        "meta": {"n": n, "w_peak": wP, "w_bg": wB, "off_px": off, "center_px": mid},
        "slices": {"peak": (rp.start, rp.stop, cp.start, cp.stop),
                   "LL":   (rll.start, rll.stop, cll.start, cll.stop)}
    }

    return mean_map, out

# --- 2) Plot using the exact same slices ---
def plot_px_windows(mean_map, slices, path, title="Mean map with pixel-scaled P2LL windows"):
    n = mean_map.shape[0]
    fig, ax = plt.subplots(figsize=(4,4))
    ax.imshow(mean_map, origin="upper", interpolation="nearest", cmap = 'YlOrBr', # vmin=0, vmax=4,
              extent=[-0.5, n-0.5, n-0.5, -0.5])  # pixel edges at half-integers

    def draw(r0, r1, c0, c1, **kw):
        ax.add_patch(Rectangle((c0 - 0.5, r0 - 0.5), (c1 - c0), (r1 - r0),
                               fill=False, **kw))

    r0p,r1p,c0p,c1p = slices["peak"]
    r0l,r1l,c0l,c1l = slices["LL"]

    draw(r0p, r1p, c0p, c1p, linewidth=2)                 # peak (solid)
    draw(r0l, r1l, c0l, c1l, linewidth=1.5, linestyle="--")  # LL (dashed)
    # ax.text(c0l - 0.3, r0l - 0.3, "LL", fontsize=8)

    ax.set_title(title)
    ax.set_xlabel("cols"); ax.set_ylabel("rows")
    plt.colorbar(ax.images[0], ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(path)
    # return fig, ax


#### The functions from cooltools that I edited!
def get_n_pixels(bad_bin_mask, window=10, ignore_diags=2):
    """
    Calculate the number of "good" pixels in a diamond at each bin.
    """

    N = len(bad_bin_mask)
    n_pixels = np.zeros(N)
    loc_bad_bin_mask = np.zeros(N, dtype=bool)
    for i_shift in range(0, window):
        for j_shift in range(0, window):
            if i_shift + j_shift < ignore_diags:
                continue

            loc_bad_bin_mask[:] = False
            if i_shift == 0:
                loc_bad_bin_mask |= bad_bin_mask
            else:
                loc_bad_bin_mask[i_shift:] |= bad_bin_mask[:-i_shift]
            if j_shift == 0:
                loc_bad_bin_mask |= bad_bin_mask
            else:
                loc_bad_bin_mask[:-j_shift] |= bad_bin_mask[j_shift:]

            n_pixels[i_shift : (-j_shift if j_shift else None)] += (
                1 - loc_bad_bin_mask[i_shift : (-j_shift if j_shift else None)]
            )
    return n_pixels


def insul_diamond(
    pixel_query,
    bins,
    window=10,
    ignore_diags=2,
    norm_by_median=True,
    clr_weight_name="weight",
):
    """
    Calculates the insulation score of a Hi-C interaction matrix.

    Parameters
    ----------
    pixel_query : RangeQuery object <TODO:update description>
        A table of Hi-C interactions. Must follow the Cooler columnar format:
        bin1_id, bin2_id, count, balanced (optional)).
    bins : pandas.DataFrame
        A table of bins, is used to determine the span of the matrix
        and the locations of bad bins.
    window : int
        The width (in bins) of the diamond window to calculate the insulation
        score.
    ignore_diags : int
        If > 0, the interactions at separations < `ignore_diags` are ignored
        when calculating the insulation score. Typically, a few first diagonals
        of the Hi-C map should be ignored due to contamination with Hi-C
        artifacts.
    norm_by_median : bool
        If True, normalize the insulation score by its NaN-median.
    clr_weight_name : str or None
        Name of balancing weight column from the cooler to use.
        Using raw unbalanced data is not supported for insulation.
    """
    lo_bin_id = bins.index.min()
    hi_bin_id = bins.index.max() + 1
    N = hi_bin_id - lo_bin_id
    sum_counts = np.zeros(N)
    sum_balanced = np.zeros(N)

    if clr_weight_name is None:
        # define n_pixels
        n_pixels = get_n_pixels(
            np.repeat(False, len(bins)), window=window, ignore_diags=ignore_diags
        )
    else:
        # calculate n_pixels
        n_pixels = get_n_pixels(
            bins[clr_weight_name].isnull().values,
            window=window,
            ignore_diags=ignore_diags,
        )
        # define transform - balanced and raw ('count') for now
        weight1 = clr_weight_name + "1"
        weight2 = clr_weight_name + "2"
        transform = lambda p: p["count"] * p[weight1] * p[weight2]

    for chunk_dict in pixel_query.read_chunked():
        chunk = pd.DataFrame(chunk_dict, columns=["bin1_id", "bin2_id", "count"])
        diag_pixels = chunk[chunk.bin2_id - chunk.bin1_id <= (window - 1) * 2]

        if clr_weight_name:
            diag_pixels = cooler.annotate(diag_pixels, bins[[clr_weight_name]])
            diag_pixels["balanced"] = transform(diag_pixels)
            valid_pixel_mask = ~diag_pixels["balanced"].isnull().values

        i = diag_pixels.bin1_id.values - lo_bin_id
        j = diag_pixels.bin2_id.values - lo_bin_id

        for i_shift in range(0, window):
            for j_shift in range(0, window):
                if i_shift + j_shift < ignore_diags:
                    continue

                mask = (
                    (i + i_shift == j - j_shift)
                    & (i + i_shift < N)
                    & (j - j_shift >= 0)
                )

                sum_counts += np.bincount(
                    i[mask] + i_shift, diag_pixels["count"].values[mask], minlength=N
                )

                if clr_weight_name:
                    sum_balanced += np.bincount(
                        i[mask & valid_pixel_mask] + i_shift,
                        diag_pixels["balanced"].values[mask & valid_pixel_mask],
                        minlength=N,
                    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        if clr_weight_name:
            score = sum_balanced / n_pixels
        else:
            score = sum_counts / n_pixels

        if norm_by_median:
            score /= np.nanmedian(score)

    return score, n_pixels, sum_balanced, sum_counts

@pool_decorator
def calculate_insulation_score(
    clr,
    window_bp,
    view_df=None,
    ignore_diags=None,
    min_dist_bad_bin=0,
    is_bad_bin_key="is_bad_bin",
    append_raw_scores=False,
    chunksize=20000000,
    clr_weight_name="weight",
    verbose=False,
    nproc=1,
    map_functor=map,
):
    """Calculate the diamond insulation scores for all bins in a cooler.

    Parameters
    ----------
    clr : cooler.Cooler
        A cooler with balanced Hi-C data.
    window_bp : int or list of integers
        The size of the sliding diamond window used to calculate the insulation
        score. If a list is provided, then a insulation score if calculated for each
        value of window_bp.
    view_df : bioframe.viewframe or None
        Viewframe for independent calculation of insulation scores for regions
    ignore_diags : int | None
        The number of diagonals to ignore. If None, equals the number of
        diagonals ignored during IC balancing.
    min_dist_bad_bin : int
        The minimal allowed distance to a bad bin to report insulation score.
        Fills bins that have a bad bin closer than this distance by nans.
    is_bad_bin_key : str
        Name of the output column to store bad bins
    append_raw_scores : bool
        If True, append columns with raw scores (sum_counts, sum_balanced, n_pixels)
        to the output table.
    clr_weight_name : str or None
        Name of the column in the bin table with weight.
        Using unbalanced data with `None` will avoid masking "bad" pixels.
    verbose : bool
        If True, report real-time progress.
    nproc : int, optional
        How many processes to use for calculation. Ignored if map_functor is passed.
    map_functor : callable, optional
        Map function to dispatch the matrix chunks to workers.
        If left unspecified, pool_decorator applies the following defaults: if nproc>1 this defaults to multiprocess.Pool;
        If nproc=1 this defaults the builtin map.

    Returns
    -------
    ins_table : pandas.DataFrame
        A table containing the insulation scores of the genomic bins
    """

    if view_df is None:
        view_df = make_cooler_view(clr)
    else:
        # Make sure view_df is a proper viewframe
        try:
            _ = is_compatible_viewframe(
                view_df,
                clr,
                check_sorting=True,
                raise_errors=True,
            )
        except Exception as e:
            raise ValueError("view_df is not a valid viewframe or incompatible") from e

    # check if cooler is balanced
    if clr_weight_name:
        try:
            _ = is_cooler_balanced(clr, clr_weight_name, raise_errors=True)
        except Exception as e:
            raise ValueError(
                f"provided cooler is not balanced or {clr_weight_name} is missing"
            ) from e

    bin_size = clr.info["bin-size"]
    # check if ignore_diags is valid
    if ignore_diags is None:
        try:
            ignore_diags = clr._load_attrs(
                clr.root.rstrip("/") + f"/bins/{clr_weight_name}"
            )["ignore_diags"]
        except:
            raise ValueError(
                f"ignore_diags not provided, and not found in cooler balancing weights {clr_weight_name}"
            )
    elif isinstance(ignore_diags, int):
        pass  # keep it as is
    else:
        raise ValueError(f"ignore_diags must be int or None, got {ignore_diags}")

    if np.isscalar(window_bp):
        window_bp = [window_bp]
    window_bp = np.array(window_bp, dtype=int)

    bad_win_sizes = window_bp % bin_size != 0
    if np.any(bad_win_sizes):
        raise ValueError(
            f"The window sizes {window_bp[bad_win_sizes]} has to be a multiple of the bin size {bin_size}"
        )

    # Calculate insulation score for each region separately.
    # Using try-clause to close mp.Pool properly
    # Apply get_region_insulation:
    job = partial(
        _get_region_insulation,
        clr,
        is_bad_bin_key,
        clr_weight_name,
        chunksize,
        window_bp,
        min_dist_bad_bin,
        ignore_diags,
        append_raw_scores,
        verbose,
    )
    ins_region_tables = map_functor(job, view_df[["chrom", "start", "end", "name"]].values)

    ins_table = pd.concat(ins_region_tables)
    return ins_table

def _get_region_insulation(
    clr,
    is_bad_bin_key,
    clr_weight_name,
    chunksize,
    window_bp,
    min_dist_bad_bin,
    ignore_diags,
    append_raw_scores,
    verbose,
    region,
):
    """
    Auxilary function to make calculate_insulation_score parallel.
    """

    # XXX -- Use a delayed query executor
    nbins = len(clr.bins())
    selector = CSRSelector(
        clr.open("r"), shape=(nbins, nbins), field="count", chunksize=chunksize
    )

    # Convert window sizes to bins:
    bin_size = clr.info["bin-size"]
    window_bins = window_bp // bin_size

    # Parse region and set up insulation table for the region:
    chrom, start, end, name = region
    region = [chrom, start, end]
    region_bins = clr.bins().fetch(region)
    ins_region = region_bins[["chrom", "start", "end"]].copy()
    ins_region.loc[:, "region"] = name
    ins_region[is_bad_bin_key] = (
        region_bins[clr_weight_name].isnull() if clr_weight_name else False
    )

    if verbose:
        logging.info(f"Processing region {name}")

    if min_dist_bad_bin:
        ins_region = ins_region.assign(
            dist_bad_bin=numutils.dist_to_mask(ins_region[is_bad_bin_key])
        )

    # XXX --- Create a delayed selection
    c0, c1 = clr.extent(region)
    region_query = selector[c0:c1, c0:c1]

    for j, win_bin in enumerate(window_bins):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            # XXX -- updated insul_diamond
            ins_track, n_pixels, sum_balanced, sum_counts = insul_diamond(
                region_query,
                region_bins,
                window=win_bin,
                ignore_diags=ignore_diags,
                clr_weight_name=clr_weight_name,
            )
            # ins_track[ins_track == 0] = np.nan
            # ins_track = np.log2(ins_track)

        # ins_track[~np.isfinite(ins_track)] = np.nan

        ins_region[f"insulation_score_{window_bp[j]}"] = ins_track
        ins_region[f"n_valid_pixels_{window_bp[j]}"] = n_pixels

        if min_dist_bad_bin:
            mask_bad = ins_region.dist_bad_bin.values < min_dist_bad_bin
            ins_region.loc[mask_bad, f"insulation_score_{window_bp[j]}"] = np.nan

        if append_raw_scores:
            ins_region[f"sum_counts_{window_bp[j]}"] = sum_counts
            ins_region[f"sum_balanced_{window_bp[j]}"] = sum_balanced

    return ins_region


def find_boundaries(
    ins_table,
    min_frac_valid_pixels=0.66,
    min_dist_bad_bin=0,
    log2_ins_key="insulation_score_{WINDOW}",
    n_valid_pixels_key="n_valid_pixels_{WINDOW}",
    is_bad_bin_key="is_bad_bin",
):
    """Call insulating boundaries.

    Find all local minima of the log2(insulation score) and calculate their
    chromosome-wide topographic prominence.

    Parameters
    ----------
    ins_table : pandas.DataFrame
        A bin table with columns containing log2(insulation score),
        annotation of regions (required),
        the number of valid pixels per diamond and (optionally) the mask
        of bad bins. Normally, this should be an output of calculate_insulation_score.
    view_df : bioframe.viewframe or None
        Viewframe for independent boundary calls for regions
    min_frac_valid_pixels : float
        The minimal fraction of valid pixels in a diamond to be used in
        boundary picking and prominence calculation.
    min_dist_bad_bin : int
        The minimal allowed distance to a bad bin to be used in boundary picking.
        Ignore bins that have a bad bin closer than this distance.
    log2_ins_key, n_valid_pixels_key : str
        The names of the columns containing log2_insulation_score and
        the number of valid pixels per diamond. When a template
        containing `{WINDOW}` is provided, the calculation is repeated
        for all pairs of columns matching the template.

    Returns
    -------
    ins_table : pandas.DataFrame
        A bin table with appended columns with boundary prominences.
    """

    if min_dist_bad_bin:
        ins_table = pd.concat(
            [
                df.assign(dist_bad_bin=numutils.dist_to_mask(df[is_bad_bin_key]))
                for region, df in ins_table.groupby("region")
            ]
        )

    if "{WINDOW}" in log2_ins_key:
        windows = set()
        for col in ins_table.columns:
            m = re.match(log2_ins_key.format(WINDOW=r"(\d+)"), col)
            if m:
                windows.add(int(m.groups()[0]))
    else:
        windows = set([None])

    min_valid_pixels = {
        win: ins_table[n_valid_pixels_key.format(WINDOW=win)].max()
        * min_frac_valid_pixels
        for win in windows
    }

    dfs = []
    index_name = ins_table.index.name  # Store the name of the index and soring order
    sorting_order = ins_table.index.values
    ins_table.index.name = "sorting_index"
    ins_table.reset_index(drop=False, inplace=True)
    for region, df in ins_table.groupby("region"):
        df = df.sort_values(["start"])  # Force sorting by the bin start coordinate
        for win in windows:
            mask = (
                df[n_valid_pixels_key.format(WINDOW=win)].values
                >= min_valid_pixels[win]
            )

            if min_dist_bad_bin:
                mask &= df.dist_bad_bin.values >= min_dist_bad_bin

            ins_track = df[log2_ins_key.format(WINDOW=win)].values[mask]
            poss, proms = peaks.find_peak_prominence(-ins_track)
            ins_prom_track = np.zeros_like(ins_track) * np.nan
            ins_prom_track[poss] = proms

            if win is not None:
                bs_key = f"boundary_strength_{win}"
            else:
                bs_key = "boundary_strength"

            df[bs_key] = np.nan
            df.loc[mask, bs_key] = ins_prom_track

        dfs.append(df)

    df = pd.concat(dfs)
    df = df.set_index("sorting_index")  # Restore original sorting order and name
    df.index.name = index_name
    df = df.loc[sorting_order, :]
    return df


def _insul_diamond_dense(mat, window=10, ignore_diags=2, norm_by_median=True):
    """
    Calculates the insulation score of a Hi-C interaction matrix.

    Parameters
    ----------
    mat : numpy.array
        A dense square matrix of Hi-C interaction frequencies.
        May contain nans, e.g. in rows/columns excluded from the analysis.

    window : int
        The width of the window to calculate the insulation score.

    ignore_diags : int
        If > 0, the interactions at separations < `ignore_diags` are ignored
        when calculating the insulation score. Typically, a few first diagonals
        of the Hi-C map should be ignored due to contamination with Hi-C
        artifacts.

    norm_by_median : bool
        If True, normalize the insulation score by its NaN-median.

    Returns
    -------
    score : ndarray
        an array with normalized insulation scores for provided matrix
    """
    if ignore_diags:
        mat = mat.copy()
        for i in range(-ignore_diags + 1, ignore_diags):
            numutils.set_diag(mat, np.nan, i)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        N = mat.shape[0]
        score = np.nan * np.ones(N)
        for i in range(0, N):
            lo = max(0, i + 1 - window)
            hi = min(i + window, N)
            # nanmean of interactions to reduce the effect of bad bins
            score[i] = np.nanmean(mat[lo : i + 1, i:hi])
        if norm_by_median:
            score /= np.nanmedian(score)
    return score



def _find_insulating_boundaries_dense(
    clr,
    window_bp=100000,
    view_df=None,
    clr_weight_name="weight",
    min_dist_bad_bin=2,
    ignore_diags=None,
):
    """Calculate the diamond insulation scores and call insulating boundaries.

    Parameters
    ----------
    clr : cooler.Cooler
        A cooler with balanced Hi-C data. Balancing weights are required
        for the detection of bad_bins.
    window_bp : int
        The size of the sliding diamond window used to calculate the insulation
        score.
    view_df : bioframe.viewframe or None
        Viewframe for independent calculation of insulation scores for regions
    clr_weight_name : str
        Name of the column in bin table that stores the balancing weights.
    min_dist_bad_bin : int
        The minimal allowed distance to a bad bin. Do not calculate insulation
        scores for bins having a bad bin closer than this distance.
    ignore_diags : int
        The number of diagonals to ignore. If None, equals the number of
        diagonals ignored during IC balancing.

    Returns
    -------
    ins_table : pandas.DataFrame
        A table containing the insulation scores of the genomic bins and
        the insulating boundary strengths.
    """

    if view_df is None:
        view_df = make_cooler_view(clr)
    else:
        # Make sure view_df is a proper viewframe
        try:
            _ = is_compatible_viewframe(
                view_df,
                clr,
                check_sorting=True,
                raise_errors=True,
            )
        except Exception as e:
            raise ValueError("view_df is not a valid viewframe or incompatible") from e

    bin_size = clr.info["bin-size"]

    # check if cooler is balanced
    if clr_weight_name:
        try:
            _ = is_cooler_balanced(clr, clr_weight_name, raise_errors=True)
        except Exception as e:
            raise ValueError(
                f"provided cooler is not balanced or {clr_weight_name} is missing"
            ) from e

    # check if ignore_diags is valid
    if ignore_diags is None:
        ignore_diags = clr._load_attrs(
            clr.root.rstrip("/") + f"/bins/{clr_weight_name}"
        )["ignore_diags"]
    elif isinstance(ignore_diags, int):
        pass  # keep it as is
    else:
        raise ValueError(f"provided ignore_diags {ignore_diags} is not int or None")

    window_bins = window_bp // bin_size

    if window_bp % bin_size != 0:
        raise ValueError(
            f"The window size ({window_bp}) has to be a multiple of the bin size {bin_size}"
        )

    ins_region_tables = []
    for chrom, start, end, name in view_df[["chrom", "start", "end", "name"]].values:
        region = [chrom, start, end]
        ins_region = clr.bins().fetch(region)[["chrom", "start", "end"]]
        is_bad_bin = np.isnan(clr.bins().fetch(region)[clr_weight_name].values)
        # extract dense Hi-C heatmap for a given "region"
        m = clr.matrix(balance=clr_weight_name).fetch(region)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            ins_track = _insul_diamond_dense(m, window_bins, ignore_diags)
            ins_track[ins_track == 0] = np.nan
            ins_track = np.log2(ins_track)

        bad_bin_neighbor = np.zeros_like(is_bad_bin)
        for i in range(0, min_dist_bad_bin):
            if i == 0:
                bad_bin_neighbor = bad_bin_neighbor | is_bad_bin
            else:
                bad_bin_neighbor = bad_bin_neighbor | np.r_[[True] * i, is_bad_bin[:-i]]
                bad_bin_neighbor = bad_bin_neighbor | np.r_[is_bad_bin[i:], [True] * i]

        ins_track[bad_bin_neighbor] = np.nan
        ins_region["bad_bin_masked"] = bad_bin_neighbor

        ins_track[~np.isfinite(ins_track)] = np.nan

        ins_region[f"insulation_score_{window_bp}"] = ins_track

        poss, proms = peaks.find_peak_prominence(-ins_track)
        ins_prom_track = np.zeros_like(ins_track) * np.nan
        ins_prom_track[poss] = proms
        ins_region[f"boundary_strength_{window_bp}"] = ins_prom_track
        ins_region[f"boundary_strength_{window_bp}"] = ins_prom_track

        ins_region_tables.append(ins_region)

    ins_table = pd.concat(ins_region_tables)
    return ins_table


def insulation(
    clr,
    window_bp,
    view_df=None,
    ignore_diags=None,
    clr_weight_name="weight",
    min_frac_valid_pixels=0.66,
    min_dist_bad_bin=0,
    threshold="Li",
    append_raw_scores=False,
    chunksize=20000000,
    verbose=False,
    nproc=1,
):
    """Find insulating boundaries in a contact map via the diamond insulation score.

    For a given cooler, this function (a) calculates the diamond insulation score track,
    (b) detects all insulating boundaries, and (c) removes weak boundaries via an automated
    thresholding algorithm.

    Parameters
    ----------
    clr : cooler.Cooler
        A cooler with balanced Hi-C data.
    window_bp : int or list of integers
        The size of the sliding diamond window used to calculate the insulation
        score. If a list is provided, then a insulation score if done for each
        value of window_bp.
    view_df : bioframe.viewframe or None
        Viewframe for independent calculation of insulation scores for regions
    ignore_diags : int | None
        The number of diagonals to ignore. If None, equals the number of
        diagonals ignored during IC balancing.
    clr_weight_name : str
        Name of the column in the bin table with weight
    min_frac_valid_pixels : float
        The minimal fraction of valid pixels in a diamond to be used in
        boundary picking and prominence calculation.
    min_dist_bad_bin : int
        The minimal allowed distance to a bad bin to report insulation score.
        Fills bins that have a bad bin closer than this distance by nans.
    threshold : "Li", "Otsu" or float
        Rule used to threshold the histogram of boundary strengths to exclude weak
        boundaries. "Li" or "Otsu" use corresponding methods from skimage.thresholding.
        Providing a float value will filter by a fixed threshold
    append_raw_scores : bool
        If True, append columns with raw scores (sum_counts, sum_balanced, n_pixels)
        to the output table.
    verbose : bool
        If True, report real-time progress.
    nproc : int, optional
        How many processes to use for calculation

    Returns
    -------
    ins_table : pandas.DataFrame
        A table containing the insulation scores of the genomic bins
    """
    # Create view:
    print('running my version of insulation!')
    if view_df is None:
        # full chromosomes:
        view_df = make_cooler_view(clr)
    else:
        # Make sure view_df is a proper viewframe
        try:
            _ = is_compatible_viewframe(
                view_df,
                clr,
                # must be sorted for pairwise regions combinations
                # to be in the upper right of the heatmap
                check_sorting=True,
                raise_errors=True,
            )
        except Exception as e:
            raise ValueError("view_df is not a valid viewframe or incompatible") from e

    if threshold == "Li":
        thresholding_func = lambda x: x >= threshold_li(x)
    elif threshold == "Otsu":
        thresholding_func = lambda x: x >= threshold_otsu(x)
    else:
        try:
            thr = float(threshold)
            thresholding_func = lambda x: x >= thr
        except ValueError:
            raise ValueError(
                "Insulating boundary strength threshold can be Li, Otsu or a float"
            )
    # Calculate insulation score:
    ins_table = calculate_insulation_score(
        clr,
        view_df=view_df,
        window_bp=window_bp,
        ignore_diags=ignore_diags,
        min_dist_bad_bin=min_dist_bad_bin,
        append_raw_scores=append_raw_scores,
        clr_weight_name=clr_weight_name,
        chunksize=chunksize,
        verbose=verbose,
        nproc=nproc,
    )

    # Find boundaries:
    ins_table = find_boundaries(
        ins_table,
        min_frac_valid_pixels=min_frac_valid_pixels,
        min_dist_bad_bin=min_dist_bad_bin,
    )
    for win in window_bp:
        strong_boundaries = thresholding_func(
            ins_table[f"boundary_strength_{win}"].values
        )
        ins_table[f"is_boundary_{win}"] = strong_boundaries

    return ins_table




if __name__ == '__main__':
    clr_path = sys.argv[1]
    cell_type = sys.argv[2] 
    genome = sys.argv[3]
    dom_path= sys.argv[4]
    dots_path = sys.argv[5]
    output_dir = sys.argv[6]
    sample_id = sys.argv[7]
    print("clr_path: " + clr_path)
    print("cell_type: " + cell_type)
    print("genome: " + genome)
    print("dom_path: " + dom_path)
    print("dots_path: " + dots_path)
    print("output_dir: " + output_dir)
    print("sample_id: " + sample_id)

    summary_rows = []
    aggregate_res=[10000, 25000, 50000, 100000] #Do the analysis for all the resolutions, so I don't need to redo it later

    for res in aggregate_res:
        print(res)
        sample_row = {"sample": sample_id, "resolution": res}
        res_dir_TAD = dom_path + '/' + str(res//1000) + 'kb/'
        ('loading clr object...')
        #### Plotting the Aggregated TADs
        clr = cooler.Cooler(clr_path + '::resolutions/' + str(res))
        # print(f'chromosomes: {clr.chromnames}, binsize: {clr.binsize}')

        # Auto-detect if cooler uses 'chr' prefix
        clr_has_chr = any(str(c).startswith('chr') for c in clr.chromnames)
        print(f"DEBUG: cooler has 'chr' prefix: {clr_has_chr}")

        if 'mm10' in genome:
            mm10_chromsizes = bioframe.fetch_chromsizes("mm10")
            mm10 = bioframe.make_viewframe(mm10_chromsizes)
            if not clr_has_chr:  # strip 'chr' only if cooler doesn't have it
                mm10.name = [x.replace('chr', '') for x in mm10.name]
                mm10.chrom = [str(x.replace('chr', '')) for x in mm10.chrom]
            mm10 = mm10[mm10.name.isin(clr.chromnames)]
            mm10["chrom"] = pd.Categorical(mm10["chrom"], categories=clr.chromnames, ordered=True)
            gene_df = mm10.sort_values("chrom").reset_index(drop=True)
            print('Calculated expected for ' + genome + '...')
        if 'hg19' in genome:
            hg19_chromsizes = bioframe.fetch_chromsizes('hg19')
            hg19_cens = bioframe.fetch_centromeres('hg19')
            hg19_arms = bioframe.make_chromarms(hg19_chromsizes, hg19_cens)
            if not clr_has_chr:  # strip 'chr' only if cooler doesn't have it
                hg19_arms.name = [x.replace('chr', '') for x in hg19_arms.name]
                hg19_arms.chrom = [str(x.replace('chr', '')) for x in hg19_arms.chrom]
            hg19_arms = hg19_arms[hg19_arms.chrom.isin(clr.chromnames)]
            hg19_arms["chrom"] = pd.Categorical(hg19_arms["chrom"], categories=clr.chromnames, ordered=True)
            gene_df = hg19_arms.sort_values(["chrom", "start"]).reset_index(drop=True)
            print('Calculated expected for ' + genome + '...')

        # DEBUG: Check chromosome naming before expected_cis
        print(f"DEBUG: clr.chromnames: {list(clr.chromnames)[:5]}")
        print(f"DEBUG: gene_df shape: {gene_df.shape}")
        print(f"DEBUG: gene_df chrom unique: {list(gene_df['chrom'].unique())[:5]}")

        expected = cooltools.expected_cis(clr, view_df=gene_df, nproc=2, chunksize=1_000_000)

        bs_stats = {}
        if res != 5000:
            print('Creating aggregate TAD plot...')
            col_names = ["chrom", "start", "end", "name", "score", "strand"]

            # Extract reference name from path (last directory name)
            # dom_path is like: /path/to/TAD/merged_ngn2_treatment_vehicle_noXY/
            tad_ref_name = os.path.basename(os.path.normpath(dom_path))

            # DEBUG: print path components
            print(f"DEBUG: dom_path = {dom_path}", flush=True)
            print(f"DEBUG: tad_ref_name = {tad_ref_name}", flush=True)
            print(f"DEBUG: res_dir_TAD = {res_dir_TAD}", flush=True)

            # Build file paths using actual naming convention
            # Files are named: TAD_{tad_ref_name}_{res}kb_min{res*3}kb_max{res*5}kb_step{res}kb_thr001_fdr001_{domains|boundaries}.bed
            dom_path_res = res_dir_TAD + 'TAD_' + tad_ref_name + '_' + str(res//1000) + 'kb_min' + str(res*3//1000) + 'kb_max' + str(res*5//1000) + 'kb_step' + str(res//1000) + 'kb_thr001_fdr001_domains.bed'
            bound_path_res = dom_path_res.replace('domains', 'boundaries')
            print(f"DEBUG: bound_path_res = {bound_path_res}", flush=True)

            try:
                df_domain = pd.read_csv(bound_path_res, sep="\t", header=None, names=col_names, comment="#", index_col = False, dtype={"chrom": str})
                if not clr_has_chr:  # strip 'chr' only if cooler doesn't have it
                    df_domain['chrom'] = df_domain['chrom'].astype(str).str.replace('^chr', '', regex=True)
            except FileNotFoundError:
                print(f"WARNING: Boundary file not found at {bound_path_res}, skipping TAD analysis...")
                df_domain = pd.DataFrame()  # empty dataframe to skip

            # Skip if boundary file is empty or not found
            if df_domain.empty:
                print(f"WARNING: No boundaries found at {res} resolution, skipping TAD analysis...")
            else:
                ("Pileup...")

                # DEBUG: Check chromosome naming
                print(f"DEBUG: clr.chromnames: {list(clr.chromnames)[:5]}")
                print(f"DEBUG: df_domain shape: {df_domain.shape}")
                print(f"DEBUG: df_domain chrom unique: {list(df_domain['chrom'].unique())[:5]}")
                print(f"DEBUG: df_domain chrom dtype: {df_domain['chrom'].dtype}")
                print(f"DEBUG: gene_df shape: {gene_df.shape}")
                print(f"DEBUG: gene_df chrom unique: {list(gene_df['chrom'].unique())[:5]}")

                pup_dom = coolpup.pileup(clr, df_domain, features_format='bed', view_df=gene_df, local=True,
                                        expected_df=expected, rescale= False, flank=res*30)
                aggregate_plot1(pup_dom['data'], 'coolwarm', res, output_dir, cell_type, 'TADs')
                plt.close('all')

                df_domain = pd.read_csv(dom_path_res, sep="\t", header=None, names=col_names, comment="#", index_col = False, dtype={"chrom": str})
                if not clr_has_chr:  # strip 'chr' only if cooler doesn't have it
                    df_domain['chrom'] = df_domain['chrom'].astype(str).str.replace('^chr', '', regex=True)
                # print('TAD sanity check plots...')
                # for i in range(1,3):
                #   a, b = random_pair(len(df_domain))
                #   print(a,b)
                #   TAD_sanity_check(clr, df_domain, df_domain.iloc[a].chrom, start = df_domain.iloc[a].start, end = df_domain.iloc[a].end + 3_000_000, resolution = res, path = output_dir)
                # print('TADs done! Starting Loop Identification!')

                print('Calculating boundary strength...')
                merged_bound = pd.read_csv(bound_path_res, sep="\t", header=None, names=["chrom","start","end","name","score","strand"])
                if not clr_has_chr:  # strip 'chr' only if cooler doesn't have it
                    merged_bound['chrom'] = merged_bound['chrom'].astype(str).str.replace('^chr', '', regex=True)
                # boundary_strength_df = get_boundary_strength(clr, merged_bound, 300_000)
                bs_stats = boundary_strength_summary(clr, merged_bound, output_dir + '/' + f"{sample_id}_{res}")
            # boundary_strength_df.to_csv(output_dir + cell_type + '_' + sample_id + '_' + res + '_boundary_strength_df.csv')
        sample_row.update(bs_stats)

        lp_stats = {}
        if res in [10000, 25000]:
            res_dir_Loops = dots_path + '/' + str(res//1000) + 'kb/'
            Loop_path_res= res_dir_Loops + str(res) + '_Loops_df.csv'
            try:
                dots_df = pd.read_csv(Loop_path_res)
                if dots_df.empty:
                    print(f"WARNING: No loops found at {res} resolution, skipping loop analysis...")
                else:
                    print("Creating aggregate Loop plot")
                    loops_bedpe = dots_df[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]
                    if not clr_has_chr:  # strip 'chr' only if cooler doesn't have it
                        loops_bedpe['chrom1'] = loops_bedpe['chrom1'].astype(str).str.replace('^chr', '', regex=True)
                        loops_bedpe['chrom2'] = loops_bedpe['chrom2'].astype(str).str.replace('^chr', '', regex=True)
                    # pup_loops = coolpup.pileup(clr, loops_bedpe, features_format='bedpe', view_df=gene_df, expected_df=expected, rescale = True)
                    pup_loops = coolpup.pileup(clr, loops_bedpe, features_format='bedpe', view_df=gene_df, expected_df=expected) #not rescaling the pixels
                    aggregate_plot2(pup_loops['data'], 'YlOrBr', res, output_dir, cell_type, 'Loops')
                    plt.close('all')
                    print('Calculating Loop stats...')
                    pup_loops.to_csv(output_dir + '/Loop_pup_df_' + str(res) + '.csv')
                    mean_map, res_out = p2ll_px(pup_loops)
                    plot_px_windows(mean_map, res_out["slices"], output_dir + '/P2LL_ratio_plot_' + str(res) + '.png')

                    #print the matrix 9x9 into a text file
                    # print('Loop sanity check plots...')
                    # for i in range(1,6):
                    #     a, b = random_pair(len(dots_df))
                    #     print(a,b)
                    #     loops_sanity_check(clr, dots_df, dots_df.iloc[a].chrom1, dots_df.iloc[a].start1, dots_df.iloc[a].start1 + 3_000_000, res, path = output_dir)
                    # print('Loops done!')
                    lp_stats = {"P2LL": res_out["P2LL"], "peak": res_out["peak"], "LL": res_out["LL"]}
            except FileNotFoundError:
                print(f"WARNING: Loop file not found at {Loop_path_res}, skipping loop analysis...")
        sample_row.update(lp_stats)
        summary_rows.append(sample_row)

    summary_df = (pd.DataFrame(summary_rows)
                .set_index(["sample", "resolution"])   # rows: sample ▸ res
                .sort_index())
    out_path = output_dir + '/' + f"{sample_id}_metrics.csv"
    print('Saving summary table to: ' + out_path)
    summary_df.to_csv(out_path, sep=",")
    print("saved →", out_path)
