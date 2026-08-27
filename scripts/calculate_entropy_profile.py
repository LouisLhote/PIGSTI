import pysam
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import sys

# --- Handle both named and positional Snakemake outputs ---
try:
    bam = snakemake.input.bam
except AttributeError:
    bam = snakemake.input[0]

try:
    plot = snakemake.output.plot
except AttributeError:
    plot = snakemake.output[0]

try:
    summary = snakemake.output.mean
except AttributeError:
    try:
        summary = snakemake.output.summary
    except AttributeError:
        summary = snakemake.output[1]

try:
    entropy_100_file = snakemake.output.entropy_100bp
    entropy_1000_file = snakemake.output.entropy_1000bp
except AttributeError:
    summary_dir = os.path.dirname(summary)
    summary_base = os.path.basename(summary).replace(".mean_entropy.txt", "")
    entropy_100_file = os.path.join(summary_dir, f"{summary_base}.entropy_100bp.txt")
    entropy_1000_file = os.path.join(summary_dir, f"{summary_base}.entropy_1000bp.txt")

# --- Optional: Check BAM index exists ---
bam_index = bam + ".bai"
if not os.path.exists(bam_index):
    sys.stderr.write(f"ERROR: Missing BAM index file: {bam_index}\n")
    sys.exit(1)

# --- Pathopipe method: Relative entropy of read start positions ---
# As described in Nature paper (Sikora et al. 2025):
# "we calculate the frequency of read alignments with their start positions 
# falling within windows along the reference assembly, using two different 
# window sizes (100 and 1,000 bp). The obtained frequency vector is converted 
# into Shannon information entropy, and normalized using the maximum entropy 
# attainable if the same total number of reads were evenly distributed across 
# the windows."

bamfile = pysam.AlignmentFile(bam, "rb")
ref_name = bamfile.references[0]
ref_len = bamfile.lengths[0]

# Count read start positions per base
starts = np.zeros(ref_len, dtype=int)
for read in bamfile.fetch(ref_name):
    if not read.is_unmapped:
        start = read.reference_start
        if start < ref_len:
            starts[start] += 1

bamfile.close()

def calculate_relative_entropy(window_size, starts, ref_len):
    """
    Calculate relative entropy using pathopipe method:
    - Partition genome into non-overlapping windows of size window_size
    - Count read starts per window
    - Compute Shannon entropy of window frequency distribution
    - Normalize by log2(number of windows)
    
    Returns: relative_entropy (0-1 scale)
    """
    # Partition genome into non-overlapping windows
    num_windows = int(np.ceil(ref_len / window_size))
    window_counts = np.zeros(num_windows, dtype=int)
    
    # Count read starts per window
    for pos in range(ref_len):
        if starts[pos] > 0:
            window_idx = pos // window_size
            if window_idx < num_windows:
                window_counts[window_idx] += starts[pos]
    
    total_reads = np.sum(window_counts)
    
    # If no reads, return 0.0
    if total_reads == 0:
        return 0.0
    
    # If only one window has reads, entropy is 0 (all reads in one bin)
    if np.sum(window_counts > 0) <= 1:
        return 0.0
    
    # Calculate probabilities: p_i = count_i / total_reads
    p = window_counts.astype(float) / total_reads
    
    # Shannon entropy: H = -Σ(p_i * log2(p_i)) for p_i > 0
    entropy = -np.sum(p[p > 0] * np.log2(p[p > 0]))
    
    # Normalization: maximum entropy if reads evenly distributed across all windows
    # H_max = log2(number of windows)
    max_entropy = math.log2(num_windows)
    
    # Relative entropy: H / H_max
    relative_ent = entropy / max_entropy if max_entropy > 0 else 0.0
    
    return relative_ent

# Calculate relative entropy for both window sizes (as in pathopipe)
rel_entropy_100 = calculate_relative_entropy(100, starts, ref_len)
rel_entropy_1000 = calculate_relative_entropy(1000, starts, ref_len)

# Write both values to summary file (for backward compatibility, also write the 1000bp value as "mean_entropy")
with open(summary, "w") as f:
    f.write(f"{rel_entropy_1000:.4f}\n")  # Keep backward compatibility with existing code
    f.write(f"# Relative entropy (100bp windows): {rel_entropy_100:.4f}\n")
    f.write(f"# Relative entropy (1000bp windows): {rel_entropy_1000:.4f}\n")

with open(entropy_100_file, "w") as f:
    f.write(f"{rel_entropy_100:.4f}\n")

with open(entropy_1000_file, "w") as f:
    f.write(f"{rel_entropy_1000:.4f}\n")

# Generate plot showing both entropy values
# Use MetBrewer style colors
colors = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']
color1 = colors[0]  # Blue for 100bp
color2 = colors[1]  # Red for 1000bp

fig, ax = plt.subplots(figsize=(10, 6))

# Create bar plot for both entropy values
metrics = ['100 bp windows', '1000 bp windows']
values = [rel_entropy_100, rel_entropy_1000]
bar_colors = [color1, color2]

bars = ax.bar(metrics, values, color=bar_colors, alpha=0.8, edgecolor='black', linewidth=1.5)

# Add threshold lines
ax.axhline(y=0.9, color='#228833', linestyle='--', linewidth=1.5, alpha=0.7, label='Bacteria threshold (0.9)')
ax.axhline(y=0.7, color='#EE6677', linestyle='--', linewidth=1.5, alpha=0.7, label='Virus threshold (0.7)')

# Add value labels on bars
for i, (bar, val) in enumerate(zip(bars, values)):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
            f'{val:.4f}', ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.set_ylim([0, 1.1])
ax.set_ylabel('Relative Entropy', fontsize=12, fontweight='bold')
ax.set_title('Relative Entropy of Read Start Positions (Pathopipe Method)', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, linestyle='--', axis='y')
ax.legend(loc='upper right', fontsize=10)
ax.set_xticklabels(metrics, fontsize=11)

plt.tight_layout()
plt.savefig(plot, dpi=300, facecolor='white')
plt.close()

sys.stderr.write(f"[INFO] Relative entropy (100bp windows): {rel_entropy_100:.4f}\n")
sys.stderr.write(f"[INFO] Relative entropy (1000bp windows): {rel_entropy_1000:.4f}\n")
