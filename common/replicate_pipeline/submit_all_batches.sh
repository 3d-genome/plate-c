#!/usr/bin/env bash
# Submit batches 1-4 sequentially, each waiting for the previous to complete
#
# Usage: bash submit_all_batches.sh
#
# Submits: 1_qc_passed_samples_other_noiPSC.csv
#          2_qc_passed_samples_other_noiPSC.csv
#          3_qc_passed_samples_other_noiPSC.csv
#          4_qc_passed_samples_other_noiPSC.csv

set -euo pipefail

PREV_JOB=""

for i in 1 2 3 4; do
    CSV_NAME="${i}_qc_passed_samples_other_noiPSC.csv"

    echo "================================================"
    echo "Submitting batch $i: $CSV_NAME"
    echo "================================================"

    if [[ -z "$PREV_JOB" ]]; then
        # First batch - no dependency
        bash submit_pipeline.sh "$CSV_NAME"
    else
        # Subsequent batches - wait for previous to complete
        echo "Waiting for job $PREV_JOB to complete..."
        bash submit_pipeline.sh "$CSV_NAME" --dependency=afterany:$PREV_JOB
    fi

    # Capture the last job ID from submit_pipeline.sh output
    # This assumes submit_pipeline.sh prints job IDs
    echo ""
done

echo "All batches submitted!"
