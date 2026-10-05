# Upstream pipeline: reads → contacts (SLURM job arrays)

These steps run once per well (Plate-C) or per cell (Easy Dip-C) on a SLURM cluster and produce the files that every attribute folder starts from. Each `*_job_array.sh` takes a list file (one sample folder or FASTQ per line) and processes line `$SLURM_ARRAY_TASK_ID`, e.g. `sbatch --array=1-384 map_mouse_job_array.sh fastq_r1_list.txt`.

| Step | Script | Output (per sample folder) |
| --- | --- | --- |
| 1. Align | `map_mouse_job_array.sh` / `map_human_job_array.sh`: `bwa mem -5SP` (BWA 0.7.17) to GRCm38 (mm10) or hs37d5 (hg19) | `aln.sam.gz` |
| 2. Contacts | `hickit_2d_bulk_remove_blacklist_{mouse,human}_job_array.sh`: `hickit.js sam2seg` → `chronly` → ENCODE blacklist filter (`bedflt`) → `hickit --dup-dist=1` | `contacts_unisex.seg.gz`, `contacts_unisex.pairs.gz`, `contacts_unisex.info` (sample, raw contacts, duplicate %, contacts, % intrachromosomal) |
| 3. Dip-C format | `dip-c_prep_job_array.sh` (`hickit_pairs_to_con.sh`) | `contacts_unisex.con.gz` |
| 4. Merge | `merge_pairs_plateC.sh header cellnames out.pairs`: merges replicates (`hickit --dup-dist=0`) | merged `.pairs.gz` |
| 5. `.hic` (optional, for viewing) | `juicer_{mouse,human}.sh`: Juicer Tools 1.22.01 `pre` + `addNorm` | `.hic` |
| 6. QC table | `extract_contacts_unisex.sh processed folders.txt out.tsv` | per-sample contacts summary and R1 read counts |

Tool paths (`/home/users/tttt/tools/...`), genome indices and the blacklist BED are set at the top of each script. Edit them for your system. Merged `.pairs.gz` files are converted to `.mcool` with `../cooler_preprocessing/pairs_to_mcool.sh`.
