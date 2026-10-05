#!/usr/bin/env bash
# Plate-C demo: runs the replicate-level analysis of the Figure 2 attributes on
# the small datasets in demo/data and writes results to demo/output.
# Run from the repository root:  bash demo/run_demo.sh
set -euo pipefail

cd "$(dirname "$0")/.."
OUT=demo/output
mkdir -p "$OUT"

echo "[1/5] Chromosome intermingling (% interchromosomal contacts vs vehicle)"
python 01_chromosome_intermingling/intermingling_vs_vehicle.py \
    demo/data/intermingling "$OUT/01_intermingling_vs_vehicle.tsv"

echo "[2/5] Contact distance distribution"
python 02_contact_distance_distribution/plot_distance_distribution.py \
    demo/data/distance_histograms/*_treatment_DMSO.log10_histogram_matrix.tsv \
    demo/data/distance_histograms/*_treatment_SGI-1027.log10_histogram_matrix.tsv \
    --colors DMSO=black SGI-1027=#dc0000 \
    --output "$OUT/02_distance_distribution.svg" \
    --summary "$OUT/02_distance_distribution_means.tsv"

echo "[3/5] Compartment strength and A-B difference vs vehicle"
python common/compare_to_vehicle.py \
    --input demo/data/saddle_strength_example.tsv \
    --treatment-col drug --vehicle DMSO \
    --value-cols saddle_strength_extent10_log2 relative_strength_extent10_log2 \
    --output "$OUT/03_04_compartment_strength_ab_difference_vs_vehicle.tsv"

echo "[4/5] Locus-level scA/B track plot"
python plot_contact_maps/plot_ab_tracks_heatmap.py \
    demo/data/demo_scab.bedgraph chr11 0 120000000 "$OUT/05_scab_track_chr11.png"

echo "[5/5] Boundary strength (aggregate insulation) and loop strength (P2LL) vs vehicle"
python common/compare_to_vehicle.py \
    --input demo/data/loop_insulation_example.tsv \
    --treatment-col treatment_index --vehicle DMSO \
    --value-cols log2_mean_aggr_ins_all P2LL \
    --output "$OUT/06_07_boundary_and_loop_strength_vs_vehicle.tsv"

echo "Done. Compare $OUT with demo/expected_output."
