# 5. Locus-level scA/B profile (1-Mb bins)

**Definition.** Single-cell A/B (scA/B) value of each 1-Mb locus: the CpG density of a locus's contact partners, computed with [`dip-c`](https://github.com/tanlongzhi/dip-c) (`dip-c color2`, CpG frequency track, 1-Mb bins) and merged across replicates with `dip-c mgcolor` into `*.cpg_b1m.color2s` matrices.

**Computing scA/B (SLURM).**
- `dip-c_mouse_job_array.sh` / `dip-c_human_job_array.sh`: `dip-c color2 -b1000000 -H -c {mm10,hg19}.cpg.1m.txt` on each `contacts_unisex.con.gz` → `cpg_b1m.color2` (dip-c runs on Python 2.7 with numpy and scipy). `run_color2_perfile.sh` does the same for a single `.con.gz` file.
- Merge samples into a matrix with `dip-c mgcolor <sample>/cpg_b1m.color2 ... > <experiment>.cpg_b1m.color2s`.
- `colors_d_to_matlab_{mouse,human}.sh` converts chromosome names to numbers for MATLAB (mouse X/Y → 20/21, human X/Y → 23/24).

**MATLAB code (R2021b or later; Statistics and Machine Learning Toolbox).**
- `analyze_platec_mouse_granule_primary.m`: Plate-C granule-cell screen. PCA/UMAP of scA/B profiles, Ward clustering, and differential scA/B per 1-Mb locus (two-sided equal-variance `ttest2`, BH FDR via `mafdr(...,'BHFDR',true)`; thresholds set by `max_fdr` and `min_diff` in the script).
- `analyze_dipc_mouse_hdaci_in_vivo.m`: the same analysis for Easy Dip-C single cells from in vivo experiments.
- `compare_delta_scab.m`: correlation and total-least-squares slope between two sets of differential scA/B profiles (`*.all_b1m_diff.txt`).
- `gene_position.mm10.vM25.midpoint_matlab_protein_coding.txt`: gene midpoints used to map loci to genes.

Input paths at the top of each script point to the authors' folders. Change `folder` and the `readcell(...)` paths to your copies of the Zenodo files.

Differential tracks can be drawn with `plot_contact_maps/plot_ab_tracks_heatmap.py`.
