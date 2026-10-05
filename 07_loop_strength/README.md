# 7. Loop strength (P2LL)

**Definition.** Peak-to-lower-left ratio, per replicate, at 25-kb resolution. The replicate's contact map is piled up (expected-normalized) over all reference loops. P2LL is the mean of the log2 pileup in the central (peak) window minus the mean in the lower-left (background) window, each window a centred quarter-width square rounded up to an odd number of pixels, the LL window in the lower-left corner (`p2ll_px`).

**Pipeline.**
1. **Reference loops (merged vehicle maps).** `Loop_finder.py` calls loops with `cooltools.dots` (≤ 10 Mb separation) and writes `{res}_Loops_df.csv` (cis-expected over chromosome arms in human, whole chromosomes in mouse). The paper uses the **25-kb** calls. It is run by `../06_boundary_strength/Pipeline2*.sh`; `Loop_finder_plot.py` adds sanity-check plots.
2. **Per-replicate P2LL.** [`../common/replicate_pipeline/make_pup_rep.py`](../common/replicate_pipeline/make_pup_rep.py), run by `process_sample_2.sh`, piles up each replicate over the reference loops at 10 and 25 kb with `coolpuppy`. It writes `P2LL`, `peak` and `LL` to `<sample>_metrics.csv`, plus a pileup plot and a window-check image.
3. **Merged-map P2LL (optional).** `Calculate_P2LL.py <mcool_dir> <loop_dir> <hg19|mm10> <out.csv>` (via `Pipeline3.sh`) computes the same quantity for the merged maps.
4. **Statistics vs vehicle.**
```bash
python common/compare_to_vehicle.py --input <exp>_resolution-25000_loop_tad.tsv \
    --treatment-col treatment_index --vehicle DMSO \
    --value-cols P2LL --output loop_strength_vs_vehicle.tsv
```
