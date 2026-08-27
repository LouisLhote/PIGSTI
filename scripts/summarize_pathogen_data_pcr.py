#!/usr/bin/env python

import os
import sys
from pathlib import Path
import argparse
import pandas as pd
import pysam
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pigsti_naming import safe_pathogen_name

parser = argparse.ArgumentParser()
parser.add_argument("--sample", required=True, help="Biological sample ID")
parser.add_argument("--pcr", required=True, help="PCR/library ID")
parser.add_argument(
    "--evalue",
    "--escore",
    dest="evalue",
    required=True,
    help="Sample-level E-value pathogen CSV",
)
parser.add_argument("--spreadsheet", required=True)
parser.add_argument("--bam_dir", required=True, help="Per-PCR pathogen_mapping directory")
parser.add_argument("--qualimap_dir", required=True)
parser.add_argument("--damage_dir", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

spreadsheet = pd.read_csv(args.spreadsheet)
spreadsheet["Krakenuniq name"] = spreadsheet["Krakenuniq name"].str.strip()
spreadsheet["Hops name"] = spreadsheet["Hops name"].str.strip()

# Escore-based detections (sample-level, same for all PCRs of that sample)
escore_df = pd.read_csv(args.evalue)
escore_df.columns = escore_df.columns.str.strip()
escore_df["taxonomy"] = escore_df["taxonomy"].str.strip()

if "taxReads" in escore_df.columns:
    reads_col = "taxReads"
elif "# of reads" in escore_df.columns:
    reads_col = "# of reads"
else:
    raise ValueError("Missing read count column in Escore file")

detected_kraken = set(escore_df["taxonomy"])

rows = []
for _, pathogen_row in spreadsheet.iterrows():
    kraken_name = pathogen_row["Krakenuniq name"]
    if kraken_name not in detected_kraken:
        continue

    pathogen_safe = safe_pathogen_name(kraken_name)

    summary = {
        "Sample": args.sample,
        "PCR": args.pcr,
        "Pathogen": kraken_name,
        "Krakenuniq_reads": 0,
        "Escore": "NA",
        "Guellil_et_al_Evalue": "NA",
        "Detected_by_Krakenuniq": False,
        "BWA_reads": "NA",
        "Coverage": "NA",
        "Evenness": "NA",
        "Read_length_mean": "NA",
        "Read_length_median": "NA",
        "Damage_5p_CtoT": "NA",
        "Relative_entropy": "NA",
        "Breadth_ratio": "NA",  # Ratio of observed to expected breadth (as in Nature paper)
        "ANI": "NA",
    }

    # Escore (sample-level, same numbers for all PCRs)
    match = escore_df[escore_df["taxonomy"] == kraken_name]
    if not match.empty:
        summary["Detected_by_Krakenuniq"] = True
        summary["Krakenuniq_reads"] = match[reads_col].values[0]
        summary["Escore"] = match["Escore"].values[0]
        # Add Guellil et al Evalue if available
        if "Guellil_et_al_Evalue" in match.columns:
            summary["Guellil_et_al_Evalue"] = match["Guellil_et_al_Evalue"].values[0]

    # Per-PCR BAM
    bam_path = os.path.join(args.bam_dir, f"{args.pcr}_{pathogen_safe}.dedup.bam")
    if os.path.exists(bam_path):
        try:
            bam = pysam.AlignmentFile(bam_path, "rb")
            summary["BWA_reads"] = sum(1 for _ in bam.fetch(until_eof=True))
            bam.close()

            bam = pysam.AlignmentFile(bam_path, "rb")
            lengths = [
                read.query_length
                for read in bam
                if not read.is_unmapped and read.query_length is not None
            ]
            bam.close()
            if lengths:
                summary["Read_length_mean"] = round(np.mean(lengths), 2)
                summary["Read_length_median"] = round(np.median(lengths), 2)
        except Exception:
            summary["BWA_reads"] = "error"
            summary["Read_length_mean"] = "error"
            summary["Read_length_median"] = "error"

    # DamageProfiler 5p_CtoT
    dfile = os.path.join(args.damage_dir, f"damageprofiler_{pathogen_safe}", "5p_freq_misincorporations.txt")
    if os.path.exists(dfile):
        try:
            ddf = pd.read_csv(dfile, sep="\t", comment="#")
            if "C>T" in ddf.columns:
                summary["Damage_5p_CtoT"] = ddf["C>T"].iloc[0]
        except Exception:
            summary["Damage_5p_CtoT"] = "error"

    # Evenness & coverage from per-PCR Qualimap
    qfile = os.path.join(args.qualimap_dir, f"qualimap_{pathogen_safe}", "genome_results.txt")
    if os.path.exists(qfile):
        with open(qfile) as f:
            for line in f:
                line = line.strip()
                if "mean coverage" in line.lower():
                    try:
                        val = line.split("=")[-1].strip()
                        summary["Coverage"] = val
                    except Exception:
                        summary["Coverage"] = "parse_error"
                elif "There is a" in line and "coverageData >= 1X" in line:
                    try:
                        percent_str = line.split("There is a")[1].split("%")[0].strip()
                        summary["Evenness"] = f"{percent_str}%"
                    except Exception:
                        summary["Evenness"] = "parse_error"

    # Breadth ratio (observed/expected), entropy, ANI – per PCR
    breadth_ratio_file = os.path.join(args.bam_dir, f"{args.pcr}_{pathogen_safe}.breadth_ratio.txt")
    if os.path.exists(breadth_ratio_file):
        try:
            val = float(open(breadth_ratio_file).readline().strip())
            summary["Breadth_ratio"] = round(val, 4)
        except Exception:
            summary["Breadth_ratio"] = "error"

    entropy_file = os.path.join(args.bam_dir, f"{args.pcr}_{pathogen_safe}.mean_entropy.txt")
    if os.path.exists(entropy_file):
        try:
            with open(entropy_file) as f:
                summary["Relative_entropy"] = float(f.readline().strip())
        except Exception:
            summary["Relative_entropy"] = "error"

    ani_file = os.path.join(args.bam_dir, f"{args.pcr}_{pathogen_safe}.ani.txt")
    if os.path.exists(ani_file):
        try:
            with open(ani_file) as f:
                for line in f:
                    if "ANI" in line:
                        summary["ANI"] = line.strip().split("≈")[-1].strip()
                        break
        except Exception:
            summary["ANI"] = "error"

    rows.append(summary)

os.makedirs(os.path.dirname(args.output), exist_ok=True)
df = pd.DataFrame(rows)
df.to_csv(args.output, index=False)

