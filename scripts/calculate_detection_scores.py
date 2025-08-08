#!/usr/bin/env python3
"""
Pathogen Detection Scoring System
Implements a 10-point scoring system for pathogen detection based on multiple criteria:
1. E-Score threshold passed (base criterion)
2. Hops detection (at least 2 in hops output)
3. Hops edit distance (3 in hops output)
4. Hops damage (4 in hops output)
5. ANI > 0.965
6. 5′ C>T deamination ≥ 0.01
7. 3′ G>A deamination ≥ 0.01
8. Breadth ratio ≥ 0.8
9. Entropy ≥ 0.9
10. K-mer rank ≤ 2
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import glob
from pathlib import Path

def safe_name(name):
    """Convert pathogen name to safe filename format"""
    return name.replace(" ", "_").replace("/", "_")

def get_sample_ref_pairs(samples, spreadsheet_df):
    """Get sample-pathogen pairs from escore results"""
    pairs = []
    for sample in samples:
        escore_path = f"results/{sample}/Escore/pathogen/{sample}_pathogen.csv"
        if os.path.exists(escore_path):
            escore_df = pd.read_csv(escore_path)
            escore_df["taxonomy"] = escore_df["taxonomy"].str.strip()
            
            # Only add pathogens from escore that are in master spreadsheet
            for pathogen in escore_df["taxonomy"]:
                if pathogen in spreadsheet_df["Krakenuniq name"].values:
                    pairs.append((sample, pathogen))
    return pairs

def calculate_detection_score(sample, pathogen, escore_data, hops_data, bwa_data, 
                            damage_data, breadth_data, entropy_data, comparison_data):
    """
    Calculate detection score based on all 10 criteria
    Returns: (total_score, detailed_scores_dict)
    """
    score = 0
    detailed_scores = {}
    
    # 0. E-Score threshold passed (base criterion)
    escore_file = f"results/{sample}/Escore/pathogen/{sample}_pathogen.csv"
    if os.path.exists(escore_file):
        escore_df = pd.read_csv(escore_file)
        # Strip leading spaces from taxonomy column
        escore_df['taxonomy'] = escore_df['taxonomy'].str.strip()
        pathogen_row = escore_df[escore_df['taxonomy'] == pathogen]
        if not pathogen_row.empty:
            # Since these are already filtered by user thresholds, they get +1
            score += 1
            detailed_scores['escore_threshold'] = 1
        else:
            detailed_scores['escore_threshold'] = 0
    else:
        detailed_scores['escore_threshold'] = 0
    
    # 1. Hops detection (at least 2 in hops output)
    hops_key = f"{sample}_{pathogen}"
    if hops_key in hops_data and hops_data[hops_key] >= 2:
        score += 1
        detailed_scores['hops_detection'] = 1
    else:
        detailed_scores['hops_detection'] = 0
    
    # 2. Hops edit distance (3 in hops output)
    if hops_key in hops_data and hops_data[hops_key] >= 3:
        score += 1
        detailed_scores['hops_edit_distance'] = 1
    else:
        detailed_scores['hops_edit_distance'] = 0
    
    # 3. Hops damage (4 in hops output)
    if hops_key in hops_data and hops_data[hops_key] >= 4:
        score += 1
        detailed_scores['hops_damage'] = 1
    else:
        detailed_scores['hops_damage'] = 0
    
    # 4. ANI > 0.965
    ani_file = f"results/{sample}/bwa_pathogen/{sample}_{safe_name(pathogen)}.ani.txt"
    if os.path.exists(ani_file):
        try:
            ani_content = open(ani_file).read().strip()
            # Extract ANI value from the file
            if "ANI ≈" in ani_content:
                ani_value = float(ani_content.split("ANI ≈ ")[1].split("%")[0])
                if ani_value > 96.5:  # Convert percentage to decimal
                    score += 1
                    detailed_scores['ani_threshold'] = 1
                else:
                    detailed_scores['ani_threshold'] = 0
            else:
                detailed_scores['ani_threshold'] = 0
        except:
            detailed_scores['ani_threshold'] = 0
    else:
        detailed_scores['ani_threshold'] = 0
    
    # 5. 5′ C>T deamination ≥ 0.01
    damage_file = f"results/{sample}/bwa_pathogen/damageprofiler_{safe_name(pathogen)}/misincorporation.txt"
    # Check if the damageprofiler directory exists (from snakemake input)
    damage_dir = f"results/{sample}/bwa_pathogen/damageprofiler_{safe_name(pathogen)}"
    if os.path.exists(damage_file) and os.path.exists(damage_dir):
        try:
            damage_df = pd.read_csv(damage_file, sep='\t')
            # Find 5′ C>T rate (adjust column names based on your damageprofiler output)
            if 'Position' in damage_df.columns and 'C>T' in damage_df.columns:
                c_to_t_row = damage_df[damage_df['Position'] == '5']
                if not c_to_t_row.empty:
                    c_to_t_rate = float(c_to_t_row['C>T'].iloc[0])
                    if c_to_t_rate >= 0.01:
                        score += 1
                        detailed_scores['c_to_t_deamination'] = 1
                    else:
                        detailed_scores['c_to_t_deamination'] = 0
                else:
                    detailed_scores['c_to_t_deamination'] = 0
            else:
                detailed_scores['c_to_t_deamination'] = 0
        except:
            detailed_scores['c_to_t_deamination'] = 0
    else:
        detailed_scores['c_to_t_deamination'] = 0
    
    # 6. 3′ G>A deamination ≥ 0.01
    if os.path.exists(damage_file) and os.path.exists(damage_dir):
        try:
            damage_df = pd.read_csv(damage_file, sep='\t')
            if 'Position' in damage_df.columns and 'G>A' in damage_df.columns:
                g_to_a_row = damage_df[damage_df['Position'] == '3']
                if not g_to_a_row.empty:
                    g_to_a_rate = float(g_to_a_row['G>A'].iloc[0])
                    if g_to_a_rate >= 0.01:
                        score += 1
                        detailed_scores['g_to_a_deamination'] = 1
                    else:
                        detailed_scores['g_to_a_deamination'] = 0
                else:
                    detailed_scores['g_to_a_deamination'] = 0
            else:
                detailed_scores['g_to_a_deamination'] = 0
        except:
            detailed_scores['g_to_a_deamination'] = 0
    else:
        detailed_scores['g_to_a_deamination'] = 0
    
    # 7. Breadth ratio ≥ 0.8
    breadth_file = f"results/{sample}/bwa_pathogen/{sample}_{safe_name(pathogen)}.breadth_ratio.txt"
    if os.path.exists(breadth_file):
        try:
            breadth_ratio = float(open(breadth_file).read().strip())
            if breadth_ratio >= 0.8:
                score += 1
                detailed_scores['breadth_ratio'] = 1
            else:
                detailed_scores['breadth_ratio'] = 0
        except:
            detailed_scores['breadth_ratio'] = 0
    else:
        detailed_scores['breadth_ratio'] = 0
    
    # 8. Entropy ≥ 0.9
    entropy_file = f"results/{sample}/bwa_pathogen/{sample}_{safe_name(pathogen)}.mean_entropy.txt"
    if os.path.exists(entropy_file):
        try:
            entropy = float(open(entropy_file).read().strip())
            if entropy >= 0.9:
                score += 1
                detailed_scores['entropy'] = 1
            else:
                detailed_scores['entropy'] = 0
        except:
            detailed_scores['entropy'] = 0
    else:
        detailed_scores['entropy'] = 0
    
    # 9. K-mer rank ≤ 2 (from comparison results)
    comparison_file = f"results/comparison/{sample}_comparison.tsv"
    if os.path.exists(comparison_file):
        try:
            comp_df = pd.read_csv(comparison_file, sep='\t')
            pathogen_row = comp_df[comp_df['pathogen'] == pathogen]
            if not pathogen_row.empty:
                # Adjust column name based on your comparison output
                kmer_rank_col = 'kmer_rank' if 'kmer_rank' in pathogen_row.columns else 'rank'
                if kmer_rank_col in pathogen_row.columns:
                    kmer_rank = pathogen_row[kmer_rank_col].iloc[0]
                    if kmer_rank <= 2:
                        score += 1
                        detailed_scores['kmer_rank'] = 1
                    else:
                        detailed_scores['kmer_rank'] = 0
                else:
                    detailed_scores['kmer_rank'] = 0
            else:
                detailed_scores['kmer_rank'] = 0
        except:
            detailed_scores['kmer_rank'] = 0
    else:
        detailed_scores['kmer_rank'] = 0
    
    return score, detailed_scores

def main():
    # Load spreadsheet data
    spreadsheet_df = pd.read_csv("config/Pathogen_spreadsheet.csv")
    spreadsheet_df.columns = spreadsheet_df.columns.str.strip()
    
    # Load SAMPLES from the Snakefile (you might need to adjust this)
    # For now, let's get samples from the escore files
    sample_dirs = glob.glob("results/*/Escore/pathogen/*_pathogen.csv")
    SAMPLES = [os.path.basename(f).replace("_pathogen.csv", "") for f in sample_dirs]
    
    # Load hops data
    hops_df = pd.read_csv(snakemake.input.hops_results, sep='\t')
    print(f"Hops data columns: {list(hops_df.columns)}")
    print(f"Hops data shape: {hops_df.shape}")
    print(f"First few rows of hops data:")
    print(hops_df.head())
    
    # Process hops data into dictionary for easy lookup
    # The format is: pathogen_name, sample1_score, sample2_score, etc.
    hops_data = {}
    
    # Get sample names from column headers (remove the .rma6 extension)
    sample_cols = [col for col in hops_df.columns if col != 'node']
    sample_names = [col.replace('_unaligned.rma6', '') for col in sample_cols]
    
    print(f"Sample names extracted: {sample_names}")
    
    for _, row in hops_df.iterrows():
        pathogen_name = row['node'].strip('"')  # Remove quotes
        
        # Find the corresponding Krakenuniq name using the Hops name from spreadsheet
        pathogen_row = spreadsheet_df[spreadsheet_df['Hops name'] == pathogen_name]
        if not pathogen_row.empty:
            krakenuniq_name = pathogen_row['Krakenuniq name'].iloc[0]
            print(f"Processing pathogen: {pathogen_name} -> {krakenuniq_name}")
            
            # For each sample, create a key and store the score
            for i, sample in enumerate(sample_names):
                score = row[sample_cols[i]]
                key = f"{sample}_{krakenuniq_name}"
                hops_data[key] = score
                print(f"  {key}: {score}")
        else:
            print(f"Warning: No matching entry found for hops pathogen: {pathogen_name}")
    
    # Calculate scores for all sample-pathogen pairs
    all_scores = []
    detailed_scores_list = []
    
    for sample, pathogen in get_sample_ref_pairs(SAMPLES, spreadsheet_df):
        score, detailed = calculate_detection_score(
            sample, pathogen, None, hops_data, None, 
            None, None, None, None
        )
        
        all_scores.append({
            'sample': sample,
            'pathogen': pathogen,
            'total_score': score,
            'max_possible_score': 10  # Total number of criteria
        })
        
        detailed['sample'] = sample
        detailed['pathogen'] = pathogen
        detailed_scores_list.append(detailed)
    
    # Create matrices
    if all_scores:
        scores_df = pd.DataFrame(all_scores)
        matrix = scores_df.pivot(index='sample', columns='pathogen', values='total_score')
        
        # Save results
        matrix.to_csv(snakemake.output.scores_matrix)
        pd.DataFrame(detailed_scores_list).to_csv(snakemake.output.detailed_scores, index=False)
        
        # Create heatmap
        plt.figure(figsize=(15, 10))
        sns.heatmap(matrix, annot=True, cmap='RdYlBu_r', 
                    cbar_kws={'label': 'Detection Score (0-10)'})
        plt.title('Pathogen Detection Scores\n(Criteria: E-Score, Hops(3), ANI, Damage(2), Breadth, Entropy, K-mer Rank)')
        plt.tight_layout()
        plt.savefig(snakemake.output.scores_heatmap, dpi=300, bbox_inches='tight')
        
        print(f"Scoring completed for {len(all_scores)} sample-pathogen pairs")
        print(f"Average score: {scores_df['total_score'].mean():.2f}/10")
        print(f"Score distribution:")
        print(scores_df['total_score'].value_counts().sort_index())
    else:
        print("No sample-pathogen pairs found for scoring")
        # Create empty files
        pd.DataFrame().to_csv(snakemake.output.scores_matrix)
        pd.DataFrame().to_csv(snakemake.output.detailed_scores)
        plt.figure(figsize=(10, 8))
        plt.text(0.5, 0.5, 'No pathogens detected', ha='center', va='center', transform=plt.gca().transAxes)
        plt.savefig(snakemake.output.scores_heatmap)

if __name__ == "__main__":
    main()
