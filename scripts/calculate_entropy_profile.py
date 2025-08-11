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
    summary = snakemake.output.summary
except AttributeError:
    summary = snakemake.output[1]

# --- Optional: Check BAM index exists ---
bam_index = bam + ".bai"
if not os.path.exists(bam_index):
    sys.stderr.write(f"ERROR: Missing BAM index file: {bam_index}\n")
    sys.exit(1)

# --- Parameters ---
window_size = 1000
step_size = 100

bamfile = pysam.AlignmentFile(bam, "rb")
ref_name = bamfile.references[0]
ref_len = bamfile.lengths[0]

# Count read start positions
starts = np.zeros(ref_len, dtype=int)
for read in bamfile.fetch(ref_name):
    if not read.is_unmapped:
        start = read.reference_start
        if start < ref_len:
            starts[start] += 1

def relative_entropy(window_counts):
    total = np.sum(window_counts)
    if total == 0:
        return 0.0
    p = window_counts / total
    entropy = -np.sum(p[p > 0] * np.log2(p[p > 0]))
    max_entropy = math.log2(len(window_counts))
    return entropy / max_entropy if max_entropy > 0 else 0.0

positions = []
entropies = []
for i in range(0, ref_len - window_size + 1, step_size):
    window = starts[i:i + window_size]
    ent = relative_entropy(window)
    entropies.append(ent)
    positions.append(i + window_size // 2)

plt.figure(figsize=(12, 4))
plt.plot(positions, entropies, lw=1.5)
plt.title("Entropy of Read Start Positions")
plt.xlabel("Genome position (bp)")
plt.ylabel("Normalized Shannon Entropy")
plt.grid(True)
plt.tight_layout()
plt.savefig(plot, dpi=300)

mean_entropy = np.mean(entropies)
with open(summary, "w") as f:
    f.write(f"{mean_entropy:.4f}\n")
