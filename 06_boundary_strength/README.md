# 6. Boundary strength (insulation)

**Definition.** Per replicate, at 50-kb resolution: log2 of the mean diamond insulation score (300-kb window, ignore_diags = 2) across reference TAD boundaries (`log2_mean_aggr_ins_all`). Lower values mean stronger insulation. Reference boundaries are called once on the merged vehicle map of each cell type, so every replicate is scored at the same positions.

**Pipeline (SLURM; cooler, cooltools, coolpuppy, HiCExplorer).**
1. **Reference boundaries (merged vehicle maps).**
   - `../common/cooler_preprocessing/pairs_to_mcool.sh`: merged `.pairs.gz` → balanced `.mcool`, canonical autosomes only.
   - `../common/cooler_preprocessing/Pipeline1.sh`: `.mcool` → HiCExplorer `.h5`.
   - `Pipeline2.sh` / `Pipeline2_bySample.sh`: `hicFindTADs` (min/max depth 3×/5× resolution, step = resolution, `--thresholdComparisons 0.01`, `--delta 0.01`, BH FDR q ≤ 0.01) → `TAD_*_boundaries.bed`. The scripts run 5–100 kb; the paper uses the **50-kb** calls. The same scripts call loops for attribute 7.
2. **Per-replicate quantification** ([`../common/replicate_pipeline`](../common/replicate_pipeline)), started with `bash submit_pipeline.sh <sample_list.csv>`:
   - `process_sample_1.sh`: replicate `contacts_unisex.pairs.gz` → 5-kb `.cool` (`cooler cload pairs`), dropping chrM/X/Y (`create_cooler.py`).
   - `process_sample_2.sh`: `cooler zoomify` to 10/25/50/100 kb, `cooler balance`, then `make_pup_rep.py`.
   - `make_pup_rep.py` (function `boundary_strength_summary`): computes the diamond insulation score at each resolution (300-kb window, ignore 2 diagonals), using a modified copy of `cooltools.insulation` included in the script. The paper uses the 50-kb values. It keeps good bins whose diamond has > 50% valid pixels, snaps the reference boundaries to bins, and reports `n_bound`, `log2_mean_aggr_ins_all` and `log2_median_aggr_ins_all` in `<sample>_metrics.csv`. It also saves a TAD pileup (`coolpup.pileup`, local mode, flank = 30 × resolution, no rescaling, normalized by the replicate's cis-expected).
3. **Visualization.** `insulation_snippet.ipynb` plots replicate insulation profiles around reference boundaries (aggregate ± 1 Mb, isolated boundaries, single loci). Edit the paths in its first cells.
4. **Statistics vs vehicle.**
```bash
python common/compare_to_vehicle.py --input <exp>_resolution-50000_loop_tad.tsv \
    --treatment-col treatment_index --vehicle DMSO \
    --value-cols log2_mean_aggr_ins_all --output boundary_strength_vs_vehicle.tsv
```
The SLURM scripts assume `PROJECT=$HOME/research/plate_c/FINAL`, with scripts in `$PROJECT/scripts/merged/` and `$PROJECT/scripts/replicates/`. Edit `PROJECT` and the `module load` lines for your cluster.
