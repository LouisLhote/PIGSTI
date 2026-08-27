#!/usr/bin/env python

"""
Create a heatmap showing detection scores for all pathogens across all samples.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from pathlib import Path

# Get input and output from Snakemake
input_excel = snakemake.input.excel
output_png = snakemake.output.png
output_pdf = snakemake.output.pdf

# Read the merged pathogen summary
try:
    try:
        df = pd.read_excel(input_excel, sheet_name="pathogen_summary")
    except ValueError:
        df = pd.read_excel(input_excel, sheet_name=0)
except Exception as e:
    print(f"Error reading {input_excel}: {e}", file=sys.stderr)
    sys.exit(1)

if df.empty:
    print("No pathogen rows in summary; writing empty heatmap placeholder.", file=sys.stderr)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(
        0.5, 0.5, "No pathogen detections to plot",
        ha="center", va="center", fontsize=16,
    )
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.savefig(output_pdf, bbox_inches="tight")
    plt.close()
    sys.exit(0)

# Check if Score column exists (format "X/Y")
if "Score" not in df.columns:
    print("Warning: Score column not found. Creating empty heatmap.", file=sys.stderr)
    # Create empty plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.5, 0.5, "No detection scores available", 
            ha="center", va="center", fontsize=16)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.savefig(output_pdf, bbox_inches="tight")
    plt.close()
    sys.exit(0)

# Normalize key labels to strings to avoid mixed-type sorting errors
# (e.g. numeric sample IDs like 50 mixed with alphanumeric IDs like T21).
for col in ("Sample", "Pathogen"):
    if col not in df.columns:
        print(f"Error: required column '{col}' not found in pathogen summary", file=sys.stderr)
        sys.exit(1)
    df[col] = df[col].fillna("").astype(str).str.strip()


# Parse Score column (format "X/Y") to extract raw and max scores
def parse_score_parts(score_str):
    """Parse score string 'X/Y' and return (X, Y) as integers"""
    if pd.isna(score_str) or score_str == "NA" or str(score_str).strip() == "":
        return 0, 0
    try:
        s = str(score_str)
        if "/" in s:
            parts = s.split("/")
            raw = int(float(parts[0]))
            max_s = int(float(parts[1]))
        else:
            raw = int(float(s))
            # If max not provided, assume 1 (will be handled as 0/1 if needed)
            max_s = 1
        return raw, max_s
    except Exception:
        return 0, 0

df[["Raw_score", "Max_score"]] = df["Score"].apply(
    lambda s: pd.Series(parse_score_parts(s))
)

# Fill missing scores with 0
df["Raw_score"] = df["Raw_score"].fillna(0)
df["Max_score"] = df["Max_score"].fillna(0)

# Normalized score (fraction of max) for consistent coloring with per-pathogen reports
def compute_norm(row):
    raw = row["Raw_score"]
    max_s = row["Max_score"]
    if max_s is None or max_s == 0:
        return 0.0
    return float(raw) / float(max_s)

df["Norm_score"] = df.apply(compute_norm, axis=1)

# Create pivot table with raw scores (integers) for color mapping
# TRANSPOSED: Samples as columns (X-axis), Pathogens as rows (Y-axis)
pivot = df.pivot_table(
    index="Pathogen",  # Rows = Pathogens
    columns="Sample",  # Columns = Samples
    values="Raw_score",  # Use raw integer scores instead of normalized fractions
    aggfunc="first"  # Use first value if duplicates exist
)

# Create pivot table with Score strings for annotation
pivot_annot = df.pivot_table(
    index="Pathogen",
    columns="Sample",
    values="Score",
    aggfunc="first"
)

# Sort pathogens and samples for better visualization
pivot = pivot.sort_index(axis=0)  # Sort pathogens (rows)
pivot = pivot.sort_index(axis=1)  # Sort samples (columns)
pivot_annot = pivot_annot.reindex_like(pivot)  # Align annotation matrix

# Find max raw score for color scale (integer values)
max_score_in_data = pivot.max().max()
if pd.isna(max_score_in_data) or max_score_in_data == 0:
    max_score_in_data = 12  # Default max (assuming max possible score is 12 with HOPS)

# Create the heatmap with dynamic sizing (transposed dimensions)
fig, ax = plt.subplots(figsize=(max(14, len(pivot.columns) * 0.8), max(10, len(pivot) * 0.5)))

# Create custom colormap: creamy white (0) to light blue to red (cold to hot)
from matplotlib.colors import LinearSegmentedColormap
# Light blue (#E3F2FD) -> Medium blue -> Yellow -> Orange -> Red
colors = ['#FFF8E7', '#E3F2FD', '#BBDEFB', '#90CAF9', '#64B5F6', '#42A5F5', '#FFE082', '#FFB74D', '#FF8A65', '#EF5350', '#E53935', '#C62828']
n_bins = 256
cmap = LinearSegmentedColormap.from_list('cold_to_hot', colors, N=n_bins)

# Create annotation matrix (show "X/Y" format)
annot_matrix = pivot_annot.fillna("").astype(str)
annot_matrix = annot_matrix.replace("nan", "")

# Create heatmap with black strokes
sns.heatmap(
    pivot,
    annot=annot_matrix,  # Show "X/Y" scores
    fmt="s",  # String format
    cmap=cmap,
    vmin=0,
    vmax=max_score_in_data,
    square=False,
    linewidths=1.5,  # Black strokes
    linecolor="black",  # Black color for cell borders
    cbar_kws={
        "label": "Pathogen detection score",
        "shrink": 0.8,
        "pad": 0.02
    },
    ax=ax,
    annot_kws={"size": 8, "weight": "bold", "color": "black"},
    cbar=True
)

# Format colorbar to show integer ticks
cbar = ax.collections[0].colorbar
if cbar is not None:
    # Get current ticks and format as integers
    ticks = cbar.get_ticks()
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{int(t)}" for t in ticks])

# Customize the plot (transposed labels)
ax.set_title("Pathogen Detection Scores Across All Samples", 
             fontsize=16, fontweight="bold", pad=20)
ax.set_xlabel("Sample", fontsize=12, fontweight="bold")
ax.set_ylabel("Pathogen", fontsize=12, fontweight="bold")

# Rotate labels for better readability
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)

# Adjust layout
plt.tight_layout()

# Save outputs
plt.savefig(output_png, dpi=300, bbox_inches="tight", facecolor="white")
plt.savefig(output_pdf, bbox_inches="tight", facecolor="white")
plt.close()

print(f"✅ Created pathogen detection score heatmap: {output_png}")

