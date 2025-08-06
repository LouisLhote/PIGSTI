import pysam
import numpy as np
import math

bam = snakemake.input.bam
out_depth = snakemake.output.depth
out_breadth = snakemake.output.breadth
out_expected = snakemake.output.expected_breadth
out_ratio = snakemake.output.ratio

bamfile = pysam.AlignmentFile(bam, "rb")
ref_name = bamfile.references[0]
ref_len = bamfile.lengths[0]

# Calculate per-base coverage
depths = np.zeros(ref_len, dtype=int)
for pileupcolumn in bamfile.pileup(ref_name):
    depths[pileupcolumn.reference_pos] = pileupcolumn.nsegments

covered = np.sum(depths >= 1)
avg_depth = np.mean(depths)
breadth = covered / ref_len

# Expected breadth using Poisson: 1 - exp(-depth)
expected_breadth = 1 - math.exp(-avg_depth)
breadth_ratio = breadth / expected_breadth if expected_breadth > 0 else 0.0

# Write results
np.savetxt(out_depth, depths, fmt="%d")
with open(out_breadth, "w") as f: f.write(f"{breadth:.4f}\n")
with open(out_expected, "w") as f: f.write(f"{expected_breadth:.4f}\n")
with open(out_ratio, "w") as f: f.write(f"{breadth_ratio:.4f}\n")
