import pysam
import numpy as np
import math

bam = snakemake.input.bam
out_depth = snakemake.output.depth
out_breadth = snakemake.output.breadth
out_expected = snakemake.output.expected_breadth
out_ratio = snakemake.output.ratio

bamfile = pysam.AlignmentFile(bam, "rb")

# Support multi-chromosome / multi-contig references by concatenating all
# reference sequences into a single pseudo-genome coordinate system.
refs = bamfile.references
lengths = bamfile.lengths

# Compute cumulative start positions for each reference
cum_starts = {}
offset = 0
for ref_name, ref_len in zip(refs, lengths):
    cum_starts[ref_name] = offset
    offset += ref_len

genome_len = offset

# Calculate per-base coverage across the full concatenated genome
depths = np.zeros(genome_len, dtype=int)
for ref_name, ref_len in zip(refs, lengths):
    for pileupcolumn in bamfile.pileup(ref_name):
        genome_pos = cum_starts[ref_name] + pileupcolumn.reference_pos
        depths[genome_pos] = pileupcolumn.nsegments

covered = np.sum(depths >= 1)
avg_depth = np.mean(depths)
breadth = covered / genome_len

# Expected breadth using Poisson: 1 - exp(-depth)
expected_breadth = 1 - math.exp(-avg_depth)
breadth_ratio = breadth / expected_breadth if expected_breadth > 0 else 0.0

# Write results
np.savetxt(out_depth, depths, fmt="%d")
with open(out_breadth, "w") as f: f.write(f"{breadth:.4f}\n")
with open(out_expected, "w") as f: f.write(f"{expected_breadth:.4f}\n")
with open(out_ratio, "w") as f: f.write(f"{breadth_ratio:.4f}\n")
