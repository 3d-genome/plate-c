# Shared code

- `compare_to_vehicle.py`: replicate-level statistics shared by attributes 1–4, 6 and 7. For each treatment it compares technical replicates with vehicle replicates using a two-sided, equal-variance *t*-test, then applies Benjamini–Hochberg correction across treatments. Output: n, mean, s.e.m., Δ vs vehicle, p, q and −log10 q.
- `cooler_preprocessing/`: `.pairs.gz` → filtered, balanced `.mcool` (`pairs_to_mcool.sh`, `create_cooler.py`), and `.mcool` → `.h5` for HiCExplorer (`Pipeline1.sh`, `submit_Pipeline1.sh`).
- `upstream_pipeline/`: SLURM scripts from FASTQ to contacts (BWA → hickit → dip-c), merging, `.hic` conversion and the per-sample QC table. See its README.
