# 6. Boundary strength (insulation)

**Definition.** Mean log2 insulation score (`cooltools.insulation`, 25-kb bins, 125-kb window) across reference TAD boundaries. Boundaries are called once on the merged vehicle map of each cell type, so every replicate is scored at the same positions. A more negative score means stronger insulation.

**Pipeline (SLURM scripts, `cooler`/`cooltools`/`HiCExplorer`).**
1. `common/cooler_preprocessing/pairs_to_mcool.sh` converts `.pairs.gz` to a multi-resolution `.mcool` (5 kb–1 Mb, balanced), keeping only canonical autosomes (`create_cooler.py`).
2. `common/cooler_preprocessing/Pipeline1.sh` converts each resolution to HiCExplorer `.h5`.
3. `Pipeline2.sh` / `Pipeline2_bySample.sh` call boundaries on merged vehicle maps (`hicFindTADs`, min/max depth 3×/5× bin size, threshold 0.01, FDR) and loops (`../07_loop_strength/Loop_finder.py`).
4. `insulation_snippet.ipynb` computes per-replicate insulation and aggregates it around the reference boundaries (mean profile ± 1 Mb, isolated boundaries, single loci).
5. The per-replicate summary `log2_mean_aggr_ins_all` (one value per replicate in the `*_loop_tad.tsv` tables) is compared with vehicle. *The script that writes these per-replicate tables (`make_pup_rep.py`) still needs to be added here.*
```bash
python common/compare_to_vehicle.py --input <exp>_resolution-25000_loop_tad.tsv \
    --treatment-col treatment_index --vehicle DMSO \
    --value-cols log2_mean_aggr_ins_all --output boundary_strength_vs_vehicle.tsv
```
The SLURM scripts expect `PROJECT=$HOME/research/plate_c/FINAL` with the scripts copied to `$PROJECT/scripts/merged/`. Edit `PROJECT` and the `module load` lines for your cluster.
