#!/usr/bin/env python3
"""
Split a pathogen-mapped BAM into "damage" vs "no-damage" subsets.

Reads are classified as "damage" if they contain at least one end deamination
pattern:
  - within the first N aligned query bases (5' end): reference 'C' -> query 'T'
  - within the last  N aligned query bases (3' end): reference 'G' -> query 'A'

This is intended to approximate the HOPS-style "ancient (damage)" vs "default
(no-damage)" split before computing edit-distance metrics.
"""

import sys
import argparse
import pysam


def classify_read(read: pysam.AlignedSegment, ref_fa: pysam.FastaFile, window: int) -> bool:
    """
    Return True if read is "damage", False if "no-damage".
    """
    if read.is_unmapped:
        return False
    if read.query_sequence is None or len(read.query_sequence) == 0:
        return False
    if read.reference_name is None:
        return False

    # Get (query_pos, ref_pos) for aligned bases where both exist.
    aligned_pairs = [
        (qpos, rpos)
        for (qpos, rpos) in read.get_aligned_pairs(matches_only=False)
        if qpos is not None and rpos is not None
    ]
    if not aligned_pairs:
        return False

    aligned_pairs.sort(key=lambda x: x[0])  # query order (5' -> 3' in query orientation)
    first = aligned_pairs[:window]
    last = aligned_pairs[-window:] if window > 0 else []

    # 5' C->T
    has_ct = False
    for qpos, rpos in first:
        qbase = read.query_sequence[qpos].upper()
        if qbase == "N":
            continue
        rbase = ref_fa.fetch(read.reference_name, rpos, rpos + 1).upper()
        if rbase == "C" and qbase == "T":
            has_ct = True
            break

    # 3' G->A
    has_ga = False
    for qpos, rpos in last:
        qbase = read.query_sequence[qpos].upper()
        if qbase == "N":
            continue
        rbase = ref_fa.fetch(read.reference_name, rpos, rpos + 1).upper()
        if rbase == "G" and qbase == "A":
            has_ga = True
            break

    return has_ct or has_ga


def main():
    # Snakemake usage
    try:
        bam_path = snakemake.input.bam
        ref_fasta = snakemake.input.ref_fasta
        damaged_out = snakemake.output.damaged_bam
        default_out = snakemake.output.default_bam
        window = int(snakemake.params.window)
    except Exception:
        parser = argparse.ArgumentParser(description="Split BAM by end deamination patterns")
        parser.add_argument("--bam", required=True, help="Input BAM (indexed)")
        parser.add_argument("--ref-fasta", required=True, help="Reference FASTA for fetching bases")
        parser.add_argument("--damaged-bam", required=True, help="Output BAM for damage reads")
        parser.add_argument("--default-bam", required=True, help="Output BAM for non-damaged (default) reads")
        parser.add_argument("--window", type=int, default=5, help="End window size (default: 5)")
        args = parser.parse_args()
        bam_path = args.bam
        ref_fasta = args.ref_fasta
        damaged_out = args.damaged_bam
        default_out = args.default_bam
        window = int(args.window)

    if window <= 0:
        raise ValueError("--window must be > 0")

    ref_fa = pysam.FastaFile(ref_fasta)
    bam_in = pysam.AlignmentFile(bam_path, "rb")

    damaged_bam = pysam.AlignmentFile(damaged_out, "wb", header=bam_in.header)
    default_bam = pysam.AlignmentFile(default_out, "wb", header=bam_in.header)

    # Single pass over all reads.
    for read in bam_in.fetch(until_eof=True):
        # We only classify mapped reads; unmapped are ignored.
        if read.is_unmapped:
            continue
        is_damaged = classify_read(read, ref_fa, window=window)
        if is_damaged:
            damaged_bam.write(read)
        else:
            default_bam.write(read)

    damaged_bam.close()
    default_bam.close()
    bam_in.close()
    ref_fa.close()


if __name__ == "__main__":
    main()

