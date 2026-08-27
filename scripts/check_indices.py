#!/usr/bin/env python3
"""
Index Checking Script for PIGSTI Pipeline
Checks for missing BWA and Bowtie2 indices and creates status files.
"""

import pandas as pd
import yaml
import os
from pathlib import Path

def check_bwa_index(fasta_path):
    """Check if BWA index exists for a FASTA file"""
    if not fasta_path or not os.path.exists(fasta_path):
        return False
    
    # BWA index files - they keep the original extension
    index_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    
    for ext in index_extensions:
        if not os.path.exists(fasta_path + ext):
            return False
    return True

def check_bowtie2_index(fasta_path):
    """Check if Bowtie2 index exists for a FASTA file"""
    if not fasta_path or not os.path.exists(fasta_path):
        return False
    
    # Bowtie2 index files
    index_extensions = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    base_path = fasta_path.replace('.fa', '').replace('.fasta', '').replace('.fna', '')
    
    for ext in index_extensions:
        if not os.path.exists(base_path + ext):
            return False
    return True

def main():
    """Main function to check all indices"""
    try:
        print("Checking reference indices...")
        
        # Create output directories
        os.makedirs('results/index_status', exist_ok=True)
        os.makedirs('logs/index_building', exist_ok=True)
        
        # Check pathogen indices
        print("Checking pathogen BWA indices...")
        spreadsheet_df = pd.read_csv('config/Pathogen_spreadsheet.csv')
        missing_pathogen_indices = []
        
        for _, row in spreadsheet_df.iterrows():
            if pd.notna(row['bwa index']) and row['bwa index'].strip():
                ref_path = row['bwa index'].strip()
                if not check_bwa_index(ref_path):
                    missing_pathogen_indices.append(ref_path)
        
        if missing_pathogen_indices:
            print(f'Missing pathogen BWA indices: {missing_pathogen_indices}')
            with open('results/workflow/index_status/missing_pathogen_indices.txt', 'w') as f:
                for idx in missing_pathogen_indices:
                    f.write(f'{idx}\n')
            # Still create the built file to indicate checking was done
            with open('results/workflow/index_status/pathogen_indices_built.txt', 'w') as f:
                f.write('Pathogen index checking completed - some indices missing\n')
        else:
            print('All pathogen BWA indices present')
            with open('results/workflow/index_status/pathogen_indices_built.txt', 'w') as f:
                f.write('All pathogen BWA indices are present\n')
        
        # Check host indices
        print("Checking host BWA indices...")
        with open('config/config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        
        missing_host_indices = []
        if 'bwa_indices' in config:
            for host, path in config['bwa_indices'].items():
                if path and path.strip() and path.strip() != '""':
                    ref_path = path.strip()
                    if not check_bwa_index(ref_path):
                        missing_host_indices.append(ref_path)
        
        if missing_host_indices:
            print(f'Missing host BWA indices: {missing_host_indices}')
            with open('results/workflow/index_status/missing_host_indices.txt', 'w') as f:
                for idx in missing_host_indices:
                    f.write(f'{idx}\n')
            # Still create the built file to indicate checking was done
            with open('results/workflow/index_status/host_indices_built.txt', 'w') as f:
                f.write('Host index checking completed - some indices missing\n')
        else:
            print('All host BWA indices present')
            with open('results/workflow/index_status/host_indices_built.txt', 'w') as f:
                f.write('All host BWA indices are present\n')
        
        # Check mtDNA indices
        print("Checking mtDNA BWA indices...")
        missing_mtdna_indices = []
        if 'mtDNA_indices' in config:
            for host, path in config['mtDNA_indices'].items():
                if path and path.strip() and path.strip() != '""':
                    ref_path = path.strip()
                    if not check_bwa_index(ref_path):
                        missing_mtdna_indices.append(ref_path)
        
        if missing_mtdna_indices:
            print(f'Missing mtDNA BWA indices: {missing_mtdna_indices}')
            with open('results/workflow/index_status/missing_mtdna_indices.txt', 'w') as f:
                for idx in missing_mtdna_indices:
                    f.write(f'{idx}\n')
            # Still create the built file to indicate checking was done
            with open('results/workflow/index_status/mtdna_indices_built.txt', 'w') as f:
                f.write('mtDNA index checking completed - some indices missing\n')
        else:
            print('All mtDNA BWA indices present')
            with open('results/workflow/index_status/mtdna_indices_built.txt', 'w') as f:
                f.write('All mtDNA BWA indices are present\n')
        
        print("Index checking completed!")
        
    except Exception as e:
        print(f"Error during index checking: {e}")
        # Create empty status files to prevent Snakemake from failing
        os.makedirs('results/index_status', exist_ok=True)
        with open('results/workflow/index_status/pathogen_indices_built.txt', 'w') as f:
            f.write('Index checking failed\n')
        with open('results/workflow/index_status/host_indices_built.txt', 'w') as f:
            f.write('Index checking failed\n')
        with open('results/workflow/index_status/mtdna_indices_built.txt', 'w') as f:
            f.write('Index checking failed\n')
        raise

if __name__ == "__main__":
    main()
