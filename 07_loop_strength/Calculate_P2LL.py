import sys
import os
import glob
import pandas as pd
import numpy as np
import cooler
import cooltools
import bioframe
from coolpuppy import coolpup
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import warnings

# Suppress warnings
warnings.filterwarnings("ignore")

# ==============================================================================
# 1. P2LL CALCULATION (Exact copy from make_pup_rep.py)
# ==============================================================================
def p2ll_px(pup_df, peak_frac=1/4, bg_frac=1/4, off_frac=1/4,
            force_odd_peak=True, force_odd_bg=False, push_ll_to_corner=True):
    mean_map = np.log2(pup_df['data'][0])
    n = mean_map.shape[0]
    mid = (n - 1) // 2
    odd = (lambda w: w if w % 2 else w + 1)

    wP  = int(round(n * peak_frac)); wP = odd(max(1, wP)) if force_odd_peak else max(1, wP)
    wB  = int(round(n * bg_frac));   wB = odd(max(1, wB)) if force_odd_bg   else max(1, wB)
    off = max(1, int(round(n * off_frac)))

    need = (wP // 2) + (wB // 2) + 1
    off = max(off, need)
    off_max = min(n - 1 - (wB // 2) - mid, mid - (wB // 2))
    off = max(1, min(off, off_max))

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

    rp, cp   = box(mid, mid, wP)
    rll, cll = box(mid + off, mid - off, wB)

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

def plot_px_windows(mean_map, slices, path, title="Mean map with pixel-scaled P2LL windows"):
    n = mean_map.shape[0]
    fig, ax = plt.subplots(figsize=(4,4))
    ax.imshow(mean_map, origin="upper", interpolation="nearest", cmap = 'YlOrBr', 
              extent=[-0.5, n-0.5, n-0.5, -0.5]) 

    def draw(r0, r1, c0, c1, **kw):
        ax.add_patch(Rectangle((c0 - 0.5, r0 - 0.5), (c1 - c0), (r1 - r0),
                               fill=False, **kw))

    r0p,r1p,c0p,c1p = slices["peak"]
    r0l,r1l,c0l,c1l = slices["LL"]
    draw(r0p, r1p, c0p, c1p, linewidth=2)
    draw(r0l, r1l, c0l, c1l, linewidth=1.5, linestyle="--")

    ax.set_title(title)
    ax.set_xlabel("cols"); ax.set_ylabel("rows")
    plt.colorbar(ax.images[0], ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()

# ==============================================================================
# 2. GENOME HANDLING (Matches make_pup_rep.py: hg19 or mm10 only)
# ==============================================================================
def get_viewframe(genome_name, clr_chromnames):
    
    clr_has_chr = any(str(c).startswith('chr') for c in clr_chromnames)
    view_df = None

    if 'mm10' in genome_name or 'mouse' in genome_name:
        # Mouse: Use whole chromosomes
        chromsizes = bioframe.fetch_chromsizes("mm10")
        view_df = bioframe.make_viewframe(chromsizes)
        if not clr_has_chr:
            view_df.name = [x.replace('chr', '') for x in view_df.name]
            view_df.chrom = [str(x.replace('chr', '')) for x in view_df.chrom]

    elif 'hg19' in genome_name:
        # Human: Use Chromosome Arms (splits at centromere)
        chromsizes = bioframe.fetch_chromsizes('hg19')
        cens = bioframe.fetch_centromeres('hg19')
        view_df = bioframe.make_chromarms(chromsizes, cens)
        if not clr_has_chr:
            view_df.name = [x.replace('chr', '') for x in view_df.name]
            view_df.chrom = [str(x.replace('chr', '')) for x in view_df.chrom]
    
    else:
        raise ValueError(f"Genome '{genome_name}' not supported. Use 'hg19' or 'mm10'.")

    # Filter to match cooler chromosomes
    view_df = view_df[view_df.chrom.isin(clr_chromnames)]
    view_df["chrom"] = pd.Categorical(view_df["chrom"], categories=clr_chromnames, ordered=True)
    view_df = view_df.sort_values(["chrom", "start"]).reset_index(drop=True)
    return view_df

# ==============================================================================
# 3. MAIN
# ==============================================================================
def main():
    if len(sys.argv) < 5:
        print("Usage: python Calculate_P2LL.py <mcool_dir> <loop_output_dir> <genome> <output_csv>")
        sys.exit(1)

    mcool_dir = sys.argv[1]
    loop_output_dir = sys.argv[2]
    genome_name = sys.argv[3]
    output_csv = sys.argv[4]

    resolutions = [10000, 25000]
    summary_data = []

    # Get sample directories
    sample_folders = [f.path for f in os.scandir(loop_output_dir) if f.is_dir()]
    print(f"Found {len(sample_folders)} sample directories.")

    for sample_path in sample_folders:
        sample_name = os.path.basename(sample_path)
        print(f"Processing {sample_name}...")

        mcool_path = os.path.join(mcool_dir, f"{sample_name}.mcool")
        if not os.path.exists(mcool_path):
            print(f"  WARNING: mcool not found for {sample_name}. Skipping.")
            continue

        for res in resolutions:
            res_kb = res // 1000
            loop_file = os.path.join(sample_path, f"{res_kb}kb", f"{res}_Loops_df.csv")
            
            if not os.path.exists(loop_file):
                continue

            try:
                dots_df = pd.read_csv(loop_file)
                if dots_df.empty:
                    print(f"  Loop file empty for {res}. Skipping.")
                    continue

                print(f"  Analyzing {res} resolution...")
                
                clr = cooler.Cooler(f"{mcool_path}::resolutions/{res}")
                
                # Fetch viewframe using the strict logic
                view_df = get_viewframe(genome_name, clr.chromnames)
                
                # Calculate Expected (needed for pileup)
                expected = cooltools.expected_cis(clr, view_df=view_df, nproc=4, chunksize=1_000_000)

                # Standardize 'chr' prefix in dots_df
                clr_has_chr = any(str(c).startswith('chr') for c in clr.chromnames)
                if not clr_has_chr:
                     dots_df['chrom1'] = dots_df['chrom1'].astype(str).str.replace('^chr', '', regex=True)
                     dots_df['chrom2'] = dots_df['chrom2'].astype(str).str.replace('^chr', '', regex=True)

                loops_bedpe = dots_df[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]

                # Pileup & P2LL
                pup_loops = coolpup.pileup(clr, loops_bedpe, features_format='bedpe', 
                                         view_df=view_df, expected_df=expected, rescale=False)

                mean_map, res_out = p2ll_px(pup_loops)
                
                # QC Plot
                qc_path = os.path.join(sample_path, f"{res_kb}kb", f"{res}_P2LL_window_check_posthoc.png")
                plot_px_windows(mean_map, res_out["slices"], qc_path, title=f"{sample_name} {res} Loop P2LL")

                row = {
                    "sample": sample_name,
                    "resolution": res,
                    "n_loops": len(dots_df),
                    "P2LL": res_out["P2LL"],
                    "peak": res_out["peak"],
                    "LL": res_out["LL"]
                }
                summary_data.append(row)
                print(f"    -> P2LL: {res_out['P2LL']:.4f}")

            except Exception as e:
                print(f"    ERROR processing {sample_name} at {res}: {e}")

    if summary_data:
        df_summary = pd.DataFrame(summary_data)
        df_summary.to_csv(output_csv, index=False)
        print(f"Done! Summary saved to {output_csv}")
    else:
        print("No data processed.")

if __name__ == "__main__":
    main()