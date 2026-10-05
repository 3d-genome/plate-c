#!/usr/bin/env python3
"""
Count the number of contacts in each sample's pairs.gz file.
Outputs a CSV with sample_id and contact_count.
"""

import pandas as pd
import gzip
import sys
import os
from pathlib import Path

def count_contacts_in_pairs(pairs_gz_path):
    """
    Count contacts (non-comment lines) in a pairs.gz file.
    """
    if not os.path.exists(pairs_gz_path):
        return None

    count = 0
    try:
        with gzip.open(pairs_gz_path, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    count += 1
        return count
    except Exception as e:
        print(f"Error reading {pairs_gz_path}: {e}", file=sys.stderr)
        return None


def main():
    if len(sys.argv) != 2:
        print("Usage: python count_contacts.py <qc_passed_samples.csv>")
        sys.exit(1)

    csv_path = sys.argv[1]

    # Read the CSV
    df = pd.read_csv(csv_path)

    print(f"Processing {len(df)} samples from {csv_path}")
    print()

    results = []

    for idx, row in df.iterrows():
        sample_id = row['sample_id']
        sample_path = row['sample_path']

        # Construct path to pairs.gz file
        pairs_gz = os.path.join(sample_path, 'contacts_unisex.pairs.gz')

        print(f"[{idx+1}/{len(df)}] {sample_id}...", end=' ', flush=True)

        # Count contacts
        count = count_contacts_in_pairs(pairs_gz)

        if count is None:
            print(f"❌ NOT FOUND or ERROR")
            status = "missing"
        else:
            print(f"✓ {count:,} contacts")
            status = "found"

        results.append({
            'sample_id': sample_id,
            'contact_count': count,
            'pairs_file': pairs_gz,
            'status': status
        })

    # Create output DataFrame
    results_df = pd.DataFrame(results)

    # Save to CSV
    output_csv = csv_path.replace('.csv', '_contact_counts.csv')
    results_df.to_csv(output_csv, index=False)

    print()
    print(f"✓ Results saved to: {output_csv}")
    print()

    # Print summary statistics
    valid_counts = results_df[results_df['status'] == 'found']['contact_count']

    print("=== SUMMARY ===")
    print(f"Total samples: {len(results_df)}")
    print(f"Found: {(results_df['status'] == 'found').sum()}")
    print(f"Missing: {(results_df['status'] == 'missing').sum()}")

    if len(valid_counts) > 0:
        print()
        print("Contact count statistics:")
        print(f"  Min:    {valid_counts.min():,}")
        print(f"  Max:    {valid_counts.max():,}")
        print(f"  Mean:   {valid_counts.mean():,.0f}")
        print(f"  Median: {valid_counts.median():,.0f}")


if __name__ == '__main__':
    main()
