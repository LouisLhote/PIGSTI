#!/usr/bin/env python
"""
Create a comprehensive Excel summary with all samples, combining:
- Host/mtDNA metrics
- Pathogen detection and metrics
- Sample-level information
"""

import pandas as pd
import sys
import os
from pathlib import Path

# Get inputs from Snakemake
try:
    host_mtdna_file = snakemake.input.host_mtdna
    pathogen_file = snakemake.input.pathogen
    output_excel = snakemake.output.excel
except AttributeError:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--host-mtdna", required=True)
    parser.add_argument("--pathogen", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    host_mtdna_file = args.host_mtdna
    pathogen_file = args.pathogen
    output_excel = args.output

# Read input files
host_mtdna = pd.read_excel(host_mtdna_file, sheet_name="Sample_level")
# Pathogen file is an Excel file, not CSV
if pathogen_file.endswith('.xlsx') or pathogen_file.endswith('.xls'):
    try:
        pathogen = pd.read_excel(pathogen_file, sheet_name="pathogen_summary")
    except ValueError:
        pathogen = pd.read_excel(pathogen_file, sheet_name=0)
else:
    # Fallback to CSV if it's actually a CSV file
    pathogen = pd.read_csv(pathogen_file)

# Normalize column names: host/mtDNA uses "bio_sample", pathogen uses "Sample"
# Create a common "Sample" column for merging
if "bio_sample" in host_mtdna.columns and "Sample" not in host_mtdna.columns:
    host_mtdna["Sample"] = host_mtdna["bio_sample"]
elif "Sample" not in host_mtdna.columns:
    # Try to find the sample column
    sample_cols = [col for col in host_mtdna.columns if "sample" in col.lower()]
    if sample_cols:
        host_mtdna["Sample"] = host_mtdna[sample_cols[0]]
    else:
        raise ValueError(f"Could not find sample column in host/mtDNA file. Available columns: {list(host_mtdna.columns)}")

# Ensure pathogen has "Sample" column
if "Sample" not in pathogen.columns:
    # Try to find the sample column
    sample_cols = [col for col in pathogen.columns if "sample" in col.lower()]
    if sample_cols:
        pathogen["Sample"] = pathogen[sample_cols[0]]
    else:
        raise ValueError(f"Could not find sample column in pathogen file. Available columns: {list(pathogen.columns)}")

# Merge on Sample column
comprehensive = host_mtdna.merge(
    pathogen,
    on="Sample",
    how="outer",
    suffixes=("_host", "_pathogen")
)

# Reorder columns for better readability
preferred_order = [
    "Sample", "species",
    "sexing_call", "sexing_female_prob", "sexing_likelihood_ratio", "sexing_status", "sexing_plot",
    "host_endogenous_pct", "host_coverage", "mtdna_coverage",
    "raw_reads", "host_q30_reads", "duplication_rate",
    "Pathogen", "Score", "Coverage", "Mapped_reads", "ANI",
    "Relative_entropy", "Relative_entropy_100bp", "Relative_entropy_1000bp",
    "Breadth_ratio", "Edit_distance_decay_quality", "Edit_distance_decay_quality_default",
    "Genus_ranking", "Damage_5p_CtoT", "Escore", "Krakenuniq_reads",
    "Criteria_passed", "Detected_by_Krakenuniq", "Detected_by_HOPS",
]

# Get all columns and order them
existing_cols = [col for col in preferred_order if col in comprehensive.columns]
remaining_cols = [col for col in comprehensive.columns if col not in existing_cols]
final_cols = existing_cols + remaining_cols

comprehensive = comprehensive[final_cols]

# Sort by Sample, then by Pathogen (if available)
if "Pathogen" in comprehensive.columns:
    comprehensive = comprehensive.sort_values(["Sample", "Pathogen"], na_position="last")
else:
    comprehensive = comprehensive.sort_values("Sample")

# Create output directory if needed
os.makedirs(os.path.dirname(output_excel), exist_ok=True)

# Write to Excel with multiple sheets
with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
    # Sheet 1: Comprehensive (all data merged)
    comprehensive.to_excel(writer, sheet_name="All_Samples", index=False)
    
    # Sheet 2: Sample Summary (one row per sample, aggregated)
    sample_summary = host_mtdna.copy()
    if "Pathogen" in pathogen.columns:
        # Add pathogen count, best score, and list of all detected pathogens per sample
        pathogen_summary = pathogen.groupby("Sample").agg({
            "Pathogen": "count",  # Count of pathogens
            "Score": lambda x: ", ".join([str(s) for s in x.dropna().head(5)]),  # Top 5 scores
            "Coverage": "max"  # Max coverage
        }).reset_index()
        pathogen_summary.columns = ["Sample", "Pathogen_Count", "Top_Scores", "Max_Coverage"]
        
        # Add comma-separated list of all detected pathogens per sample
        pathogens_list = pathogen.groupby("Sample")["Pathogen"].apply(
            lambda x: ", ".join(sorted(set([str(p) for p in x.dropna()])))
        ).reset_index()
        pathogens_list.columns = ["Sample", "Pathogens_detected"]
        
        # Merge pathogen summary with pathogens list
        pathogen_summary = pathogen_summary.merge(pathogens_list, on="Sample", how="left")
        
        sample_summary = sample_summary.merge(pathogen_summary, on="Sample", how="left")
    
    sample_summary.to_excel(writer, sheet_name="Sample_Summary", index=False)
    
    # Sheet 3: Pathogen Summary (one row per sample-pathogen)
    if "Pathogen" in pathogen.columns:
        pathogen.to_excel(writer, sheet_name="Pathogen_Details", index=False)
    
    # Sheet 4: Host/mtDNA Summary
    host_mtdna.to_excel(writer, sheet_name="Host_mtDNA", index=False)

    # Sheet 5: Sexing (PCR-level detail + plot paths), if present
    try:
        sexing_detail = pd.read_excel(host_mtdna_file, sheet_name="Sexing")
        sexing_detail.to_excel(writer, sheet_name="Sexing", index=False)
    except ValueError:
        pass

print(f"✅ Created comprehensive summary: {output_excel}")
print(f"   - All_Samples: {len(comprehensive)} rows")
print(f"   - Sample_Summary: {len(sample_summary)} rows")
if "Pathogen" in pathogen.columns:
    print(f"   - Pathogen_Details: {len(pathogen)} rows")
print(f"   - Host_mtDNA: {len(host_mtdna)} rows")

