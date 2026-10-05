#!/usr/bin/env python3
"""
Extract all QC-passed sample paths from nested experiment folders.
Reads .txt files to get individual sample paths.
Outputs a CSV with sample information.
"""

import os
import re
from pathlib import Path
import pandas as pd

def parse_experiment_folder(folder_name):
    """
    Parse experiment folder name to extract metadata.
    Format: experiment-##_plate-c_CELLTYPE_PANEL_TIMEPOINT
    Examples:
    - experiment-01_plate-c_human_hek293_drug-panel-1a_24h
    - experiment-03_plate-c_mouse_primary-cerebellar-granule_drug-panel-1a_72h
    """
    pattern = r'experiment-(\d+)_plate-c_([^_]+)_(.+?)_drug-panel-([^_]+)_(.+)'
    match = re.match(pattern, folder_name)
    
    if match:
        exp_num, organism, cell_type, panel, timepoint = match.groups()
        return {
            'experiment': f"experiment-{exp_num}",
            'organism': organism,
            'cell_type': cell_type,
            'panel': f"drug-panel-{panel}",
            'timepoint': timepoint
        }
    return None

def read_sample_paths(txt_file):
    """
    Read a .txt file and extract all sample paths.
    Each line contains a full path to a sample directory.
    """
    samples = []
    with open(txt_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:  # Skip empty lines
                samples.append(line)
    return samples

def extract_sample_name(sample_path):
    """
    Extract the sample name from the full path.
    Example: /oak/.../parasar_251125b_plateC_16drug_eNeuron_2h_sample_025
    Returns: parasar_251125b_plateC_16drug_eNeuron_2h_sample_025
    """
    return Path(sample_path).name

def extract_sample_paths(base_dir):
    """
    Traverse the directory structure and extract all QC-passed samples.
    
    Structure:
    base_dir/
      experiment-01_plate-c_.../
        experiment-01_plate-c_..._treatment_DRUG.txt
        ...
    """
    base_path = Path(base_dir)
    samples = []
    
    # Iterate through experiment folders
    for exp_folder in sorted(base_path.iterdir()):
        if not exp_folder.is_dir():
            continue
        
        if not exp_folder.name.startswith('experiment-'):
            continue
        
        print(f"Processing: {exp_folder.name}")
        
        # Parse experiment metadata
        exp_metadata = parse_experiment_folder(exp_folder.name)
        
        # Find all .txt files matching the pattern: experiment-**_plate-c_**_treatment_**.txt
        pattern = f"{exp_folder.name}_treatment_*.txt"
        txt_files = list(exp_folder.glob(pattern))
        
        for txt_file in txt_files:
            # Extract treatment name from filename
            # Pattern: experiment-##_plate-c_..._treatment_TREATMENT-NAME.txt
            # Example: experiment-06_plate-c_human_ngn2-esc-derived_drug-panel-4_2h_treatment_4-AP.txt
            prefix = exp_folder.name + '_treatment_'
            if txt_file.stem.startswith(prefix):
                treatment = txt_file.stem.replace(prefix, '')
            else:
                # Fallback if pattern doesn't match
                treatment = txt_file.stem.replace(exp_folder.name + '_', '')
            
            # Read the sample paths from the .txt file
            sample_paths = read_sample_paths(txt_file)
            
            print(f"  - {treatment}: {len(sample_paths)} samples")
            
            # Create an entry for each sample
            for sample_path in sample_paths:
                sample_name = extract_sample_name(sample_path)
                
                sample_info = {
                    'sample_name': sample_name,
                    'sample_path': sample_path,
                    'qc_file': txt_file.name,
                    'qc_file_path': str(txt_file),
                    'treatment': treatment,
                }
                
                # Add experiment metadata if parsed successfully
                if exp_metadata:
                    sample_info.update(exp_metadata)
                else:
                    sample_info.update({
                        'experiment': exp_folder.name,
                        'organism': 'unknown',
                        'cell_type': 'unknown',
                        'panel': 'unknown',
                        'timepoint': 'unknown'
                    })
                
                samples.append(sample_info)
    
    return samples

def main():
    base_dir = "/oak/stanford/groups/tttt/users/bibudha/research/plateC/aux_data/"
    
    print(f"Scanning directory: {base_dir}")
    print("-" * 60)
    
    # Extract all sample paths
    samples = extract_sample_paths(base_dir)
    
    # Convert to DataFrame
    df = pd.DataFrame(samples)
    
    # Reorder columns for clarity
    column_order = [
        'experiment', 'organism', 'cell_type', 'panel', 'timepoint',
        'treatment', 'sample_name', 'sample_path', 'qc_file', 'qc_file_path'
    ]
    df = df[column_order]
    
    # Sort by experiment, treatment, and sample
    df = df.sort_values(['experiment', 'treatment', 'sample_name']).reset_index(drop=True)
    
    # Save to CSV
    output_file = 'qc_passed_samples.csv'
    df.to_csv(output_file, index=False)
    
    print("-" * 60)
    print(f"\n✓ Found {len(samples)} QC-passed samples")
    print(f"✓ Saved to: {output_file}")
    print(f"\nBreakdown by experiment:")
    print(df.groupby('experiment').size())
    print(f"\nBreakdown by cell type:")
    print(df.groupby('cell_type').size())
    print(f"\nBreakdown by treatment:")
    print(df.groupby('treatment').size().head(20))

if __name__ == '__main__':
    main()