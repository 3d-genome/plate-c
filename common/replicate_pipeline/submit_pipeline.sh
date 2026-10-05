#!/bin/bash
# Pipeline submission script for Hi-C processing
# Submits conversion (process_sample_1.sh) and analysis (process_sample_2.sh) jobs with dependencies
#
# Usage: bash submit_pipeline.sh <csv_filename>
# Example: bash submit_pipeline.sh qc_passed_samples_ngn2.csv

# Check if CSV filename provided
if [ $# -eq 0 ]; then
    echo "ERROR: Please provide CSV filename"
    echo "Usage: bash submit_pipeline.sh <csv_filename>"
    echo "Example: bash submit_pipeline.sh qc_passed_samples_ngn2.csv"
    exit 1
fi

CSV_NAME=$1

# Make sure logs directory exists
mkdir -p logs

# Count samples in CSV (excluding header)
CSV="$HOME/research/plate_c/FINAL/scripts/replicates/lists/$CSV_NAME"

# Check if CSV exists
if [ ! -f "$CSV" ]; then
    echo "ERROR: CSV file not found: $CSV"
    exit 1
fi

NUM_SAMPLES=$(
  awk 'NR==1{next} {sub(/\r$/,"")} $0 !~ /^[[:space:]]*$/ {n++} END{print n+0}' "$CSV"
)
LAST_INDEX=$((NUM_SAMPLES - 1))

echo "================================================"
echo "Hi-C Processing Pipeline Submission"
echo "================================================"
echo "CSV file: $CSV"
echo "Total samples: $NUM_SAMPLES"
echo "Array indices: 0-$LAST_INDEX"
echo ""

# Check if there are any samples to process
if [ $NUM_SAMPLES -le 0 ]; then
    echo "ERROR: No data rows found in CSV file (only header present)"
    echo "Please add sample data to the CSV file"
    exit 1
fi

# Submit conversion job array
echo "[1] Submitting conversion jobs (pairs → cool)..."
JOB1=$(sbatch --parsable --array=0-$LAST_INDEX --export=CSV_NAME="$CSV_NAME" process_sample_1.sh)
echo "    Job ID: $JOB1"
echo ""

# Submit analysis job array with dependency
echo "[2] Submitting analysis jobs (mcool + balance + pileup)..."
echo "    Dependency: aftercorr:$JOB1 (starts as each conversion completes)"
JOB2=$(sbatch --parsable --dependency=aftercorr:$JOB1 --array=0-$LAST_INDEX --export=CSV_NAME="$CSV_NAME" process_sample_2.sh)
echo "    Job ID: $JOB2"
echo ""

echo "================================================"
echo "Submission complete!"
echo "================================================"
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/cool_${JOB1}_0.out"
echo "  tail -f logs/analyze_${JOB2}_0.out"
echo ""
echo "Check job efficiency after completion:"
echo "  seff ${JOB1}"
echo "  seff ${JOB2}"
echo "================================================"
