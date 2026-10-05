# 7. Loop strength (P2LL)

**Definition.** Peak-to-lower-left ratio: in the log2 observed/expected pileup of all reference loops, the mean of the central (peak) window minus the mean of the lower-left (background) window. Each window is ¼ of the pileup width (`p2ll_px` in `Calculate_P2LL.py`).

**Pipeline.**
1. `Loop_finder.py` calls loops on merged vehicle maps with `cooltools.dots` (chromosome-arm view, ≤10 Mb separation) at 5–100 kb and writes `{res}_Loops_df.csv`. It is run by `../06_boundary_strength/Pipeline2*.sh`. `Loop_finder_plot.py` adds sanity-check plots.
2. `Calculate_P2LL.py <mcool_dir> <loop_dir> <hg19|mm10> <out.csv>` (run through `Pipeline3.sh`) piles up each map over its loops with `coolpuppy` (expected-normalized, no rescaling) at 10 and 25 kb, and reports P2LL, peak and LL. It also saves a QC image of the windows.
3. Per-replicate P2LL (the `P2LL` column of the `*_loop_tad.tsv` tables, computed with the same `p2ll_px` function on each replicate pileup) is compared with vehicle. *The per-replicate driver script (`make_pup_rep.py`) still needs to be added here.*
```bash
python common/compare_to_vehicle.py --input <exp>_resolution-25000_loop_tad.tsv \
    --treatment-col treatment_index --vehicle DMSO \
    --value-cols P2LL --output loop_strength_vs_vehicle.tsv
```
