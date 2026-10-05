# 2. Contact distance distribution

**Definition.** Per-replicate distribution of intrachromosomal (cis) contact distances on autosomes, for distances ≥ 1 kb. Distances are binned in log10(bp) with bin width 0.05 from 10^2.95 to 10^8.35 bp, and each replicate's histogram is normalized to sum to 1.

**Pipeline.**
1. `analyze_histogram_distance_cis_log10_job_array.sh <folder_list>` (SLURM array) writes `contacts_unisex.distance_log10_cis_histogram.txt` (bin, fraction) from each `contacts_unisex.pairs.gz`.
2. Combine the replicates of one treatment into a matrix (column 1 = bin `x`, one column per replicate):
   - `concat_treatment_log10_histogram.sh aux_data/<experiment>/<..._treatment_X>.txt` → `*.log10_histogram_matrix.tsv`;
   - `make_cluster_concatenate_log10_histogram_matrix.sh <cluster_file>.txt` builds the histograms and the matrix in one step for a sample list.
3. Test each treatment against vehicle:
   - `run_ks_permutation.sh aux_data/<experiment>`: group-level Kolmogorov–Smirnov test on ECDFs of each treatment vs vehicle. Replicate histograms are subsampled to a fixed number of contacts per replicate and ECDFs are computed on a common grid. The p-value is empirical, from permuting replicate labels between treatment and vehicle, then BH-corrected across treatments (output columns `D_obs`, `p_empirical`, `q_BH`; parameters `K_PER_REP=200000`, `GRID_N=2000`, `N_PERM=20000`, `SEED=1`, `MIN_REP=2`). As committed, the script loops over `experiment*_cluster_*.txt` files. To test treatments, switch to the commented `*_treatment_*.txt` line. `run_bh_fdr_groupKS.sh` re-runs only the BH step.
   - `run_median_distance_vs_vehicle.sh aux_data/<experiment>`: shift of the median contact distance vs vehicle (output columns `delta_median_log10`, `delta_median_bp`, `ratio_median_bp`, `p_mwu`, `q_BH`).
4. Plot replicate and mean distributions:
   ```bash
   python 02_contact_distance_distribution/plot_distance_distribution.py \
       <vehicle>.log10_histogram_matrix.tsv <treatment>.log10_histogram_matrix.tsv \
       --colors DMSO=black <drug>=#dc0000 --output distance.svg --summary means.tsv
   ```

> **To add:** the Python helpers that `run_ks_permutation.sh` and `run_median_distance_vs_vehicle.sh` call (`group_ks_ecdf_perm_with_sign.py`, `summarize_groupKS.py`, `bh_fdr_groupKS.py`, `group_median_distance.py`) are not yet in this folder.
