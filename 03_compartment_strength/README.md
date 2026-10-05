# 3. Compartment strength

**Definition.** log2[(AA + BB) / (AB + BA)] from a saddle analysis of each replicate at 1-Mb resolution.

**Pipeline.**
1. `prelim_compartments.py` balances each replicate's cooler (`cooler.balance_cooler`, mad_max = 5, min_nnz = 10) and computes GC content per bin (`bioframe.frac_gc`).
2. `compartments_and_saddlepoints.ipynb`
   - computes the compartment eigenvector (E1) of the merged vehicle map with `cooltools.eigs_cis`, GC-phased;
   - builds, for every replicate, a saddle (`cooltools.saddle`) of observed/expected contacts over 48 E1 quantile groups (2nd–98th percentile, plus the two outlier groups), using the merged E1 so that all replicates share the same A/B assignment;
   - averages the corners of the saddle (extent 10 groups) into BB, AA, AB and BA, and writes `AABBABBAvalues_*.tsv` (one row per replicate).
3. Compute `saddle_strength_extent10_log2 = log2((AA + BB) / (AB + BA))` and compare it with vehicle:
```bash
python common/compare_to_vehicle.py --input AABBABBAvalues_<exp>.tsv \
    --treatment-col drug --vehicle DMSO \
    --value-cols saddle_strength_extent10_log2 --output compartment_strength_vs_vehicle.tsv
```
