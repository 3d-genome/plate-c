# 1. Chromosome intermingling

**Definition.** Percentage of a replicate's contacts that are interchromosomal (100 − % intrachromosomal).

**Pipeline.**
1. `../common/upstream_pipeline/hickit_2d_bulk_remove_blacklist_*_job_array.sh` writes `contacts_unisex.info` for each sample. Field 5 is the % of contacts with both ends on the same chromosome.
2. Collect per-replicate values for one treatment:
   - `concat_treatment_contactloss_values.sh aux_data/<experiment>/<..._treatment_X>.txt` (Plate-C wells, via the experiment's `rename_*.txt` map) → `*.contactloss.values.txt` (one value per replicate, = 100 − field 5);
   - `make_cluster_concatenate_inter_percentage.sh <cluster_file>.txt` does the same for a list of samples (e.g. Easy Dip-C clusters) → `*.inter_percentage.txt`.
3. Compare each treatment with vehicle:
   ```bash
   python 01_chromosome_intermingling/intermingling_vs_vehicle.py <base_dir> <output.tsv>
   ```
   `<base_dir>` holds one folder per experiment, each with `*_treatment_<name>.txt` files (one sample directory per line). Exactly one file per experiment must contain `vehicle` in its name. The output gives, per treatment, mean % interchromosomal, n, a two-sided equal-variance *t*-test vs vehicle and BH FDR (per experiment).
