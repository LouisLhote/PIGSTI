#!/usr/bin/env python3
"""
Pathogen Detection Scoring System
Implements a 9-point scoring system for pathogen detection based on multiple criteria:
1. E-Score threshold passed (base criterion)
2. Hops detection (at least 2 in HOPS)
3. Hops edit distance (HOPS ≥ 3)
4. Hops damage (HOPS ≥ 4)
5. ANI > 0.965
6. 5′ C>T deamination ≥ 0.01
7. 3′ G>A deamination ≥ 0.01
8. Breadth ratio ≥ 0.8
9. Entropy ≥ 0.7
Note: The previous k-mer rank criterion was removed because rank data are not present in the comparison output.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import glob
import sys
from pathlib import Path
import re

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pigsti_naming import safe_pathogen_name as safe_name

def normalize_name(name: str) -> str:
    """Normalize pathogen names across tools for robust matching.

    - Lowercase
    - Replace underscores and hyphens with spaces
    - Remove quotes
    - Remove non-alphanumeric chars except spaces
    - Collapse multiple spaces
    """
    if name is None:
        return ""
    s = str(name).lower()
    s = s.replace("_", " ").replace("-", " ").replace("\"", "").strip()
    s = re.sub(r"[^a-z0-9\s]", "", s)
    s = re.sub(r"\s+", " ", s)
    return s

def try_read_float_file(candidate_paths):
    """Try a list of file paths and parse the first existing one as a single float.
    Returns (value or None, path_used or None)."""
    for p in candidate_paths:
        if os.path.exists(p):
            try:
                with open(p, 'r') as f:
                    txt = f.read().strip()
                    # Handle cases where the file may contain a label like 'value: X'
                    m = re.search(r"(-?\d+(?:\.\d+)?)", txt)
                    if m:
                        return float(m.group(1)), p
            except Exception:
                pass
    return None, None

def parse_damageprofiler_rates(damage_dir: str):
    """Parse damageprofiler outputs robustly to get 5' C>T and 3' G>A rates.
    Tries common filenames and column conventions.
    Returns (c_to_t_5p, g_to_a_3p). Values are floats or None.
    """
    # First, prefer explicit frequency files if present
    five_p_ct_files = [
        os.path.join(damage_dir, '5pCtoT_freq.txt'),
        os.path.join(damage_dir, '5p_freq_misincorporations.txt'),
    ]
    three_p_ga_files = [
        os.path.join(damage_dir, '3pGtoA_freq.txt'),
        os.path.join(damage_dir, '3p_freq_misincorporations.txt'),
    ]
    def read_freq_table(path, value_col_hints):
        import pandas as pd
        if not os.path.exists(path):
            return None
        try:
            df = pd.read_csv(path, sep='\t')
        except Exception:
            try:
                df = pd.read_csv(path, sep=',')
            except Exception:
                return None
        # normalize column names (remove spaces and non-alphanum)
        original_cols = list(df.columns)
        norm_map = {}
        for c in original_cols:
            norm = re.sub(r'[^A-Za-z0-9]', '', str(c)).upper()
            norm_map[norm] = c
        # find position column
        pos_col = None
        for cand in ['POS','POSITION','#POSITION']:
            if cand in norm_map:
                pos_col = norm_map[cand]
                break
        # find value column (normalize hints similarly)
        val_col = None
        hint_norms = [re.sub(r'[^A-Za-z0-9]', '', h).upper() for h in value_col_hints]
        # Try exact match first
        for hint in hint_norms:
            if hint in norm_map:
                val_col = norm_map[hint]
                break
        # Fallback: fuzzy contains (e.g., '3PG>A' -> '3PGA')
        if val_col is None:
            for norm, orig in norm_map.items():
                if any(h in norm for h in hint_norms):
                    val_col = orig
                    break
        if pos_col is None or val_col is None or df.empty:
            return None
        # Choose the first position (1) if available; else min position row
        try:
            df_pos = pd.to_numeric(df[pos_col], errors='coerce')
            if (df_pos == 1).any():
                row = df[df_pos == 1].iloc[0]
            else:
                row = df.iloc[df_pos.idxmin()]
            return float(row[val_col])
        except Exception:
            return None

    c_to_t_from_freq = None
    for fp in five_p_ct_files:
        c_to_t_from_freq = read_freq_table(fp, ['5PCT', 'CT', 'C>T', 'FREQ'])
        if c_to_t_from_freq is not None:
            break
    g_to_a_from_freq = None
    for fp in three_p_ga_files:
        g_to_a_from_freq = read_freq_table(fp, ['3PGA', 'GA', 'G>A', 'FREQ'])
        if g_to_a_from_freq is not None:
            break
    if c_to_t_from_freq is not None or g_to_a_from_freq is not None:
        return c_to_t_from_freq, g_to_a_from_freq

    # Fallback to misincorporation summary tables
    candidates = [
        os.path.join(damage_dir, 'misincorporation.txt'),
        os.path.join(damage_dir, 'misincorporation_frequency.txt'),
        os.path.join(damage_dir, 'misincorporation.tsv'),
    ]
    import pandas as pd  # local import to avoid global if script is imported elsewhere
    for fp in candidates:
        if os.path.exists(fp):
            try:
                df = pd.read_csv(fp, sep='\t')
            except Exception:
                try:
                    df = pd.read_csv(fp, sep=',')
                except Exception:
                    continue
            # Normalize column names
            df.columns = [str(c).strip() for c in df.columns]
            # Identify position column
            pos_col = None
            for cand in ['Position', 'Pos', 'position', 'pos', '#Position']:
                if cand in df.columns:
                    pos_col = cand
                    break
            # Identify deamination columns (case-insensitively)
            def find_col(substr):
                for c in df.columns:
                    c_norm = str(c).replace(' ', '').replace('-', '').upper()
                    if substr in c_norm:
                        return c
                return None
            c2t_col = find_col('C>T') or find_col('CT')
            g2a_col = find_col('G>A') or find_col('GA')
            c_to_t_5p = None
            g_to_a_3p = None
            if pos_col and (c2t_col or g2a_col):
                # Try to pick 5' as the minimal position and 3' as the maximal
                # Fall back to first/last rows if positions missing
                try:
                    positions = df[pos_col]
                    # Convert to numeric where possible
                    pos_numeric = pd.to_numeric(positions, errors='coerce')
                    if pos_numeric.notna().any():
                        idx_5 = pos_numeric.idxmin()
                        idx_3 = pos_numeric.idxmax()
                    else:
                        idx_5 = df.index.min()
                        idx_3 = df.index.max()
                except Exception:
                    idx_5 = df.index.min()
                    idx_3 = df.index.max()
                if c2t_col and idx_5 is not None:
                    try:
                        c_to_t_5p = float(df.loc[idx_5, c2t_col])
                    except Exception:
                        pass
                if g2a_col and idx_3 is not None:
                    try:
                        g_to_a_3p = float(df.loc[idx_3, g2a_col])
                    except Exception:
                        pass
            return c_to_t_5p, g_to_a_3p
    return None, None

def get_sample_ref_pairs(samples, spreadsheet_df):
    """Get sample-pathogen pairs from escore results"""
    # Use a set to avoid duplicate (sample, pathogen) pairs,
    # which break pandas pivot().
    pairs = set()
    # Precompute normalized spreadsheet Krakenuniq names for membership checks
    spreadsheet_df = spreadsheet_df.copy()
    spreadsheet_df["_krakenuniq_norm"] = spreadsheet_df["Krakenuniq name"].apply(normalize_name)
    for sample in samples:
        escore_path = f"results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv"
        if os.path.exists(escore_path):
            escore_df = pd.read_csv(escore_path)
            escore_df["taxonomy"] = escore_df["taxonomy"].astype(str).str.strip()
            escore_df["_tax_norm"] = escore_df["taxonomy"].apply(normalize_name)
            
            # Only add pathogens from escore that are in master spreadsheet, using normalized names
            allowed = set(spreadsheet_df["_krakenuniq_norm"].tolist())
            for _, row in escore_df.iterrows():
                if row["_tax_norm"] in allowed:
                    # Use the canonical Krakenuniq name from spreadsheet to avoid variant spellings
                    canonical = spreadsheet_df.loc[
                        spreadsheet_df["_krakenuniq_norm"] == row["_tax_norm"], "Krakenuniq name"
                    ].iloc[0]
                    pairs.add((sample, canonical))
    # Return a stable ordering for reproducibility
    return sorted(pairs, key=lambda x: (x[0], x[1]))

def calculate_detection_score(sample, pathogen, escore_data, hops_data, bwa_data, 
                            damage_data, breadth_data, entropy_data, comparison_data):
    """
    Calculate detection score based on all 10 criteria
    Returns: (total_score, detailed_scores_dict)
    """
    score = 0
    detailed_scores = {}
    
    # 0. E-Score threshold passed (base criterion)
    escore_file = f"results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv"
    if os.path.exists(escore_file):
        escore_df = pd.read_csv(escore_file)
        # Strip leading spaces from taxonomy column
        escore_df['taxonomy'] = escore_df['taxonomy'].astype(str).str.strip()
        escore_df['_tax_norm'] = escore_df['taxonomy'].apply(normalize_name)
        pathogen_row = escore_df[escore_df['_tax_norm'] == normalize_name(pathogen)]
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
    ani_file = f"results/pathogen/{sample}/pathogen_mapping/{sample}_{safe_name(pathogen)}.ani.txt"
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
    
    # 5. 5′ C>T and 3′ G>A deamination (ancient DNA signals)
    damage_dir = f"results/pathogen/{sample}/pathogen_mapping/damageprofiler_{safe_name(pathogen)}"
    c_to_t_rate, g_to_a_rate = parse_damageprofiler_rates(damage_dir)
    # Thresholds
    c_to_t_threshold = 0.01
    g_to_a_threshold = 0.01
    if c_to_t_rate is not None and c_to_t_rate >= c_to_t_threshold:
        score += 1
        detailed_scores['c_to_t_deamination'] = 1
    else:
        detailed_scores['c_to_t_deamination'] = 0
    if g_to_a_rate is not None and g_to_a_rate >= g_to_a_threshold:
        score += 1
        detailed_scores['g_to_a_deamination'] = 1
    else:
        detailed_scores['g_to_a_deamination'] = 0
    
    # 7. Breadth ratio ≥ 0.8
    breadth_file = f"results/pathogen/{sample}/pathogen_mapping/{sample}_{safe_name(pathogen)}.breadth_ratio.txt"
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
    # Try multiple candidate entropy files
    entropy_candidates = [
        f"results/pathogen/{sample}/pathogen_mapping/{sample}_{safe_name(pathogen)}.mean_entropy.txt",
        f"results/pathogen/{sample}/pathogen_mapping/{sample}_{safe_name(pathogen)}.entropy_mean.txt",
        f"results/pathogen/{sample}/pathogen_mapping/{sample}_{safe_name(pathogen)}.entropy.txt",
    ]
    entropy, used_entropy = try_read_float_file(entropy_candidates)
    # Slightly relaxed threshold to 0.7 to avoid over-stringency on shallow data
    entropy_threshold = 0.7
    if entropy is not None and entropy >= entropy_threshold:
        score += 1
        detailed_scores['entropy'] = 1
    else:
        detailed_scores['entropy'] = 0
    
    # 9. K-mer rank criterion removed (no rank available in comparison output)
    
    return score, detailed_scores

def main():
    # Load spreadsheet data
    spreadsheet_df = pd.read_csv("config/Pathogen_spreadsheet.csv")
    spreadsheet_df.columns = spreadsheet_df.columns.str.strip()
    # Normalize matching keys for HOPS and Krakenuniq names
    spreadsheet_df['_hops_norm'] = spreadsheet_df['Hops name'].apply(normalize_name)
    spreadsheet_df['_krakenuniq_norm'] = spreadsheet_df['Krakenuniq name'].apply(normalize_name)
    
    # Load SAMPLES from the Snakefile (you might need to adjust this)
    # For now, let's get samples from the escore files
    sample_dirs = glob.glob("results/*/evalue/pathogen/*_pathogen.csv")
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
        pathogen_name_raw = str(row['node']).strip('"')
        pathogen_name_norm = normalize_name(pathogen_name_raw)
        # Find the corresponding Krakenuniq name using normalized HOPS name from spreadsheet
        match = spreadsheet_df[spreadsheet_df['_hops_norm'] == pathogen_name_norm]
        if not match.empty:
            krakenuniq_name = match['Krakenuniq name'].iloc[0]
            print(f"Processing pathogen: {pathogen_name_raw} -> {krakenuniq_name}")
            # For each sample, create a key and store the score
            for i, sample in enumerate(sample_names):
                score = row[sample_cols[i]]
                key = f"{sample}_{krakenuniq_name}"
                hops_data[key] = score
                print(f"  {key}: {score}")
        else:
            print(f"Warning: No matching entry found for hops pathogen: {pathogen_name_raw} (normalized: {pathogen_name_norm})")
    
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
            'max_possible_score': 9  # Total number of criteria (k-mer rank removed)
        })
        
        detailed['sample'] = sample
        detailed['pathogen'] = pathogen
        detailed_scores_list.append(detailed)
    
    # Create matrices
    if all_scores:
        scores_df = pd.DataFrame(all_scores)
        # Be robust to any remaining duplicates by aggregating.
        # (Pivot fails when duplicates exist.)
        matrix = scores_df.pivot_table(
            index="sample",
            columns="pathogen",
            values="total_score",
            aggfunc="max",
        )
        
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
