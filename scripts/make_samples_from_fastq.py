#!/usr/bin/env python
"""
Generate a PIGSTI-compatible samples.tsv from a directory of FASTQ files.

Assumes filenames like:
  MYRZ180-P-UPFLEX1-NUDG1-10-189-PCR1_S56_L004_R1_001.fastq.gz
  MYRZ180-P-UPFLEX1-NUDG1-10-189-PCR1_S56_L004_R2_001.fastq.gz

Output columns:
  sample, pcr, r1, r2, RGLB, sequencing_run
"""

import argparse
import os
import re
import sys
from typing import Dict, List

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Create PIGSTI samples.tsv from FASTQ directory")
    p.add_argument(
        "--fastq-dir",
        required=True,
        help="Directory containing FASTQ(.gz) files",
    )
    p.add_argument(
        "--run-name",
        required=True,
        help='Sequencing run label to put in the "sequencing_run" column (e.g. "Novaseq39")',
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output TSV path (e.g. config/samples_novaseq39.tsv)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    fastq_dir = os.path.abspath(args.fastq_dir)
    if not os.path.isdir(fastq_dir):
        sys.stderr.write(f"ERROR: fastq-dir does not exist or is not a directory: {fastq_dir}\n")
        sys.exit(1)

    run_name = args.run_name

    # Accept multiple naming patterns:
    # 1) Illumina-style with lane: <LIBID>_Sxx_Lyyy_R[12]_001.fastq(.gz)
    #    e.g. MYRZ180-P-UPFLEX1-NUDG1-10-189-PCR1_S56_L004_R1_001.fastq.gz
    # 2) Illumina-style without lane: <LIBID>_Sxx_R[12]_001.fastq(.gz)
    #    e.g. LL033-A-MBEDTA-UDG1-15-50-PCR1_S26_R1_001.fastq.gz
    # 3) Simple R suffix with _001: <LIBID>_R[12]_001.fastq(.gz)
    #    e.g. LL034-A-MBEX1-UDG1-19-1-PCR2_R1_001.fastq.gz
    # 4) Simple R suffix: <LIBID>_R[12].fastq(.gz)
    # 5) Simple numeric suffix with _001: <LIBID>_[12]_001.fastq(.gz)
    # 6) Simple numeric suffix: <LIBID>_[12].fastq(.gz)
    pattern_illumina = re.compile(
        r"^(?P<lib>.+)_S\d+_L\d+_R(?P<read>[12])_001\.fastq(\.gz)?$"
    )
    pattern_illumina_nolane = re.compile(
        r"^(?P<lib>.+)_S\d+_R(?P<read>[12])_001\.fastq(\.gz)?$"
    )
    pattern_r_001 = re.compile(
        r"^(?P<lib>.+)_R(?P<read>[12])_001\.fastq(\.gz)?$"
    )
    pattern_r = re.compile(
        r"^(?P<lib>.+)_R(?P<read>[12])\.fastq(\.gz)?$"
    )
    pattern_num_001 = re.compile(
        r"^(?P<lib>.+)_(?P<read>[12])_001\.fastq(\.gz)?$"
    )
    pattern_num = re.compile(
        r"^(?P<lib>.+)_(?P<read>[12])\.fastq(\.gz)?$"
    )

    pairs: Dict[str, Dict[str, str]] = {}

    for fn in sorted(os.listdir(fastq_dir)):
        m = (
            pattern_illumina.match(fn)
            or pattern_illumina_nolane.match(fn)
            or pattern_r_001.match(fn)
            or pattern_r.match(fn)
            or pattern_num_001.match(fn)
            or pattern_num.match(fn)
        )
        if not m:
            continue
        lib = m.group("lib")
        read = m.group("read")
        full_path = os.path.join(fastq_dir, fn)
        pairs.setdefault(lib, {})[read] = full_path

    if not pairs:
        sys.stderr.write(
            f"WARNING: No FASTQ files matching expected pattern found in {fastq_dir}\n"
        )

    rows: List[Dict[str, str]] = []
    for lib, reads in sorted(pairs.items()):
        r1 = reads.get("1", "")
        r2 = reads.get("2", "")
        if not r1 and not r2:
            continue
        pcr_id = lib  # full library ID
        # Derive sample ID:
        # - If the library name contains '-', use the part before the first '-' (e.g. "MYRZ180-..." -> "MYRZ180")
        # - Otherwise, for Macrogen-style names like "BA27_31_1" or "KD056_40_50", use the token
        #   before the first '_' (e.g. "BA27_31_1" -> "BA27").
        if "-" in lib:
            sample_id = lib.split("-")[0]
        elif "_" in lib and re.match(r"^[A-Za-z]+[0-9]+_", lib):
            sample_id = lib.split("_")[0]
        else:
            sample_id = lib
        rows.append(
            {
                "sample": sample_id,
                "pcr": pcr_id,
                "r1": r1,
                "r2": r2,
                "RGLB": pcr_id,
                "sequencing_run": run_name,
                "source": "LOCAL",
            }
        )

    if not rows:
        sys.stderr.write("WARNING: No R1/R2 pairs detected; nothing to write.\n")

    df = pd.DataFrame(
        rows, columns=["sample", "pcr", "r1", "r2", "RGLB", "sequencing_run", "source"]
    )
    out_path = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    sys.stderr.write(f"Wrote {len(df)} rows to {out_path}\n")


if __name__ == "__main__":
    main()


