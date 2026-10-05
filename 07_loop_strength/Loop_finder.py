import warnings
warnings.filterwarnings("ignore")
from itertools import combinations
import os
import sys
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import gridspec
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
from matplotlib.ticker import EngFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import pandas as pd
from multiprocessing import Pool
import bioframe
import cooler
import cooltools
from cooltools.lib.numutils import fill_diag

from coolpuppy import coolpup
from coolpuppy import plotpup
import random
from packaging import version
bp_formatter = EngFormatter('b')

if version.parse(cooltools.__version__) < version.parse('0.5.2'):
    raise AssertionError("tutorial relies on cooltools version 0.5.2 or higher,"+
                         "please check your cooltools version and update to the latest")



def _safe_region(clr, chrom, start, end):
    """
    Clip region to chromosome bounds.
    Return (chrom, start, end) or None if invalid.
    """
    chrom = str(chrom)
    if chrom not in clr.chromsizes:
        return None

    chrom_len = int(clr.chromsizes[chrom])

    try:
        start = int(start)
        end = int(end)
    except Exception:
        return None

    start = max(0, start)
    end = min(chrom_len, end)

    if start >= end:
        return None

    return chrom, start, end


def aggregate_plot1(data, cmap, resolution, output_dir, title1, title2):

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
        cbar.set_label("log₂(mean O/E)")
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
        save_path = output_dir+'/' + title1 + '_' + str(resolution) + '_aggrplot' + '.pdf'
        save_path = save_path.replace('TADs', title2)
        plt.savefig(save_path)

def aggregate_plot2(data, cmap, resolution, output_dir, title1, title2):

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
                       # vmin=-.75, vmax=.75
                       )

        # Colorbar underneath
        cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar.set_label("log₂(mean O/E)")
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
        save_path = output_dir+'/' + title1 + '_' + str(resolution) + '_aggrplot' + '.pdf'
        save_path = save_path.replace('TAD', title2)
        plt.savefig(save_path)


def TAD_sanity_check(data, domains, chr, start, end, resolution, path):
  resolution = resolution
  region = (chr, start, end)
  start = int(start)
  end = int(end)
  if end < start:
      region = (chr, end, start)

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

  offset  = region[1]
  max_pos = (region[2] - region[1]) // resolution
  contact_matrix = np.full((max_pos, max_pos), np.nan)

  for _, row in df_region.iterrows():
      # Renamed variables to avoid overwriting the 'start'/'end' args used for filename
      bin_start = (row['start'] - offset) // resolution
      bin_end = (row['end'] - offset) // resolution

      if bin_end - bin_start > 2:
          contact_matrix[bin_start:bin_end, bin_start:bin_end] = 1
          contact_matrix[bin_start+1:bin_end-1, bin_start+1:bin_end-1] = np.nan
      else:
          contact_matrix[bin_start:bin_end, bin_start:bin_end] = 1

  im2 = pcolormesh_45deg(ax, contact_matrix, start=region[1], resolution=resolution, cmap='Blues', vmax=1, vmin=-1, alpha=0.6)
  plt.tight_layout()
  plt.savefig(path + '/' + f'{resolution}_{chr}_{start}_{end}' + '_random_TAD_sanitycheck_plot.pdf')




def loops_sanity_check(data, loops, chr, start, end, resolution, path):
  # define a region to look into as an example
  resolution = resolution
  region = (chr, start, end)
  start = int(start)
  end = int(end)
  if end < start:
      region = (chr, end, start)
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


if __name__ == '__main__':
    clr_path = sys.argv[1]
    output_dir = sys.argv[2]
    df_dom_path= sys.argv[3]
    resolution = int(sys.argv[4])
    file_name = sys.argv[5]

    print('TADs Identified!')
    print('Creating aggregate plot...')
    #### Plotting the Aggregated TADs
    clr = cooler.Cooler(clr_path + '::resolutions/' + str(resolution))
    print(f'chromosomes: {clr.chromnames}, binsize: {clr.binsize}')

    # Auto-detect if cooler uses 'chr' prefix
    clr_has_chr = any(str(c).startswith('chr') for c in clr.chromnames)
    print(f"DEBUG: cooler has 'chr' prefix: {clr_has_chr}")

    if 'mouse' in file_name:
        print('mouse')
        mm10_chromsizes = bioframe.fetch_chromsizes("mm10")
        mm10 = bioframe.make_viewframe(mm10_chromsizes)
        if not clr_has_chr:  # strip 'chr' only if cooler doesn't have it
            mm10.name = [x.replace('chr', '') for x in mm10.name]
            mm10.chrom = [str(x.replace('chr', '')) for x in mm10.chrom]
        mm10 = mm10[mm10.name.isin(clr.chromnames)]
        mm10["chrom"] = pd.Categorical(mm10["chrom"], categories=clr.chromnames, ordered=True)
        chrm_arms = mm10.sort_values("chrom").reset_index(drop=True)
        chrm_arms["chrom"] = chrm_arms["chrom"].astype(str)



    else:
        hg19_chromsizes = bioframe.fetch_chromsizes('hg19')
        hg19_cens = bioframe.fetch_centromeres('hg19')
        hg19_arms = bioframe.make_chromarms(hg19_chromsizes, hg19_cens)
        if not clr_has_chr:  # strip 'chr' only if cooler doesn't have it
            hg19_arms.name = [x.replace('chr', '') for x in hg19_arms.name]
            hg19_arms.chrom = [str(x.replace('chr', '')) for x in hg19_arms.chrom]
        hg19_arms = hg19_arms[hg19_arms.chrom.isin(clr.chromnames)]
        hg19_arms["chrom"] = pd.Categorical(hg19_arms["chrom"], categories=clr.chromnames, ordered=True)
        chrm_arms = hg19_arms.sort_values(["chrom", "start"]).reset_index(drop=True)
        chrm_arms["chrom"] = chrm_arms["chrom"].astype(str)

    print(f"DEBUG: chrm_arms shape: {chrm_arms.shape}")
    print(f"DEBUG: chrm_arms chrom unique: {list(chrm_arms['chrom'].unique())[:5]}")

    expected = cooltools.expected_cis(clr, view_df=chrm_arms, nproc=2, chunksize=1_000_000)


    if resolution in [25000, 50000, 100000]:
        col_names = ["chrom", "start", "end", "name", "score", "strand"]

        # Load both files
        df_boundaries = pd.read_csv(df_dom_path + '_boundaries.bed', sep="\t", header=None, names=col_names, comment="#", index_col=False, dtype={"chrom": str})
        df_domains = pd.read_csv(df_dom_path + '_domains.bed', sep="\t", header=None, names=col_names, comment="#", index_col=False, dtype={"chrom": str})

        # CHECK: If domains are empty, skip EVERYTHING (Pileup + Sanity Checks)
        if df_domains.empty:
             print(f"WARNING: No domains found for {resolution} resolution. Skipping TAD analysis.")

        else:
            #### 1. Plotting the Aggregated Pileup
            print("Pileup...")
            pup_dom = coolpup.pileup(clr, df_boundaries, features_format='bed', view_df=chrm_arms, local=True,
                                    expected_df=expected, rescale= False, flank=resolution*30)

            aggregate_plot1(pup_dom['data'], 'coolwarm', resolution, output_dir, file_name, 'TADs')

            #### 2. TAD sanity check plots
            print('TAD sanity check plots...')
            fixed_chrom = 'chr12' if clr_has_chr else '12'
            fixed_start = 66660000
            fixed_end   = 68050000

            try:
                TAD_sanity_check(clr, df_domains, fixed_chrom, fixed_start, fixed_end, resolution, path=output_dir)
            except Exception as e:
                print(f"Skipping fixed TAD sanity check: {e}")

            for i in range(1, 6):
                a, b = random_pair(len(df_domains))
                print(a, b)
                try:
                    TAD_sanity_check(
                        clr,
                        df_domains,
                        df_domains.iloc[a].chrom,
                        start=df_domains.iloc[a].start,
                        end=df_domains.iloc[b].end,
                        resolution=resolution,
                        path=output_dir
                    )

                except Exception as e:
                    print(f"Skipping failed TAD sanity check ({a}, {b}): {e}")
                    continue

        print('TADs done! Starting Loop Identification!')


    if resolution in [5000, 10000, 25000]:
        # Identifying Loops, aka dots
        output_dir = output_dir.replace('TAD', 'Loop')
        dots_df = cooltools.dots(clr, expected=expected, view_df=chrm_arms, max_loci_separation=10_000_000, nproc=4)
        dots_df.to_csv(output_dir + '/' + str(resolution) + '_Loops_df.csv')

        # CHECK: If no loops found, skip aggregation and sanity checks
        if dots_df.empty:
            print(f"WARNING: No loops found for {resolution} resolution. Skipping Loop analysis.")

        else:
            print("Loops Identified!")
            print("Creating aggregate plot")
            loops_bedpe = dots_df[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]

            # Pass to coolpup
            pup_loops = coolpup.pileup(clr, loops_bedpe, features_format='bedpe', view_df=chrm_arms, expected_df=expected)
            aggregate_plot2(pup_loops['data'], 'YlOrBr', resolution, output_dir, file_name, 'Loop')

            print('Loop sanity check plots...')
            fixed_chrom = 'chr12' if clr_has_chr else '12'
            fixed_start = 66660000
            fixed_end   = 68050000

            try:
                loops_sanity_check(clr, dots_df, fixed_chrom, fixed_start, fixed_end, resolution, path=output_dir)
            except Exception as e:
                print(f"Skipping fixed Loop sanity check: {e}")

            for i in range(1, 6):
                a, b = random_pair(len(dots_df))
                print(a,b)
                try:
                    loops_sanity_check(clr, dots_df, dots_df.iloc[a].chrom1, dots_df.iloc[a].start1, dots_df.iloc[b].end2, resolution, path = output_dir)
                except Exception as e:
                    print(f"Skipping failed Loop sanity check ({a}, {b}): {e}")
                    continue

        print('Loops done!')
