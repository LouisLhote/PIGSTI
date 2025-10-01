#!/usr/bin/env python

import os
import pandas as pd
import argparse
import pysam
import numpy as np

# -------------------- Argument parsing --------------------
parser = argparse.ArgumentParser()
parser.add_argument("--sample", required=True)
parser.add_argument("--escore", required=True)
parser.add_argument("--hops", required=True)
parser.add_argument("--spreadsheet", required=True)
parser.add_argument("--bam_dir", required=True)
parser.add_argument("--qualimap_dir", required=True)
parser.add_argument("--damage_dir", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

# -------------------- Load reference spreadsheet --------------------
spreadsheet = pd.read_csv(args.spreadsheet)
spreadsheet["Krakenuniq name"] = spreadsheet["Krakenuniq name"].str.strip()
spreadsheet["Hops name"] = spreadsheet["Hops name"].str.strip()

# -------------------- Load escore (Krakenuniq) --------------------
escore_df = pd.read_csv(args.escore)
escore_df.columns = escore_df.columns.str.strip()
escore_df["taxonomy"] = escore_df["taxonomy"].str.strip()

if "taxReads" in escore_df.columns:
    reads_col = "taxReads"
elif "# of reads" in escore_df.columns:
    reads_col = "# of reads"
else:
    raise ValueError("Missing read count column in Escore file")

# -------------------- Load HOPS --------------------
hops_df = pd.read_csv(args.hops, sep="\t")
hops_df.columns = hops_df.columns.str.replace('"', '').str.strip()
hops_df.rename(columns={"node": "Species"}, inplace=True)
hops_df["Species"] = hops_df["Species"].str.replace('"', '').str.strip()
sample_col = f"{args.sample}_unaligned.rma6"

# -------------------- Summarize detected pathogens only --------------------
rows = []

detected_kraken = set(escore_df["taxonomy"])
detected_hops = set()
if sample_col in hops_df.columns:
    detected_hops = set(hops_df[hops_df[sample_col] > 1]["Species"])

observed = set()
for _, row in spreadsheet.iterrows():
    if row["Krakenuniq name"] in detected_kraken:
        observed.add(row["Krakenuniq name"])
    if row["Hops name"].replace(" ", "_") in detected_hops:
        observed.add(row["Krakenuniq name"])

for _, pathogen_row in spreadsheet.iterrows():
    kraken_name = pathogen_row["Krakenuniq name"]
    if kraken_name not in observed:
        continue

    hops_name = pathogen_row["Hops name"]
    pathogen_safe = kraken_name.replace(" ", "_")
    hops_name_safe = hops_name.replace(" ", "_")

    summary = {
        "Sample": args.sample,
        "Pathogen": kraken_name,
        "Krakenuniq_reads": 0,
        "Escore": "NA",
        "Detected_by_Krakenuniq": False,
        "HOPS_score": 0,
        "Detected_by_HOPS": False,
        "BWA_reads": "NA",
        "Coverage": "NA",
        "Evenness": "NA",
        "Read_length_mean": "NA",
        "Read_length_median": "NA",
        "Damage_5p_CtoT": "NA",
        "Expected_breadth": "NA",
        "Relative_entropy": "NA",
        "ANI": "NA",
        "Entropy_plot": "NA"
    }

    # Krakenuniq match
    match = escore_df[escore_df["taxonomy"] == kraken_name]
    if not match.empty:
        summary["Detected_by_Krakenuniq"] = True
        summary["Krakenuniq_reads"] = match[reads_col].values[0]
        summary["Escore"] = match["Escore"].values[0]

    # HOPS match
    match = hops_df[hops_df["Species"] == hops_name_safe]
    if not match.empty and sample_col in match.columns:
        hops_count = match[sample_col].values[0]
        summary["HOPS_score"] = hops_count
        if hops_count > 1:
            summary["Detected_by_HOPS"] = True

    # BAM path
    bam_path = os.path.join(args.bam_dir, f"{args.sample}_{pathogen_safe}_F4_q30_sort.bam")
    if os.path.exists(bam_path):
        try:
            bam = pysam.AlignmentFile(bam_path, "rb")
            summary["BWA_reads"] = sum(1 for _ in bam.fetch(until_eof=True))
            bam.close()

            bam = pysam.AlignmentFile(bam_path, "rb")
            lengths = [read.query_length for read in bam if not read.is_unmapped and read.query_length is not None]
            bam.close()

            if lengths:
                summary["Read_length_mean"] = round(np.mean(lengths), 2)
                summary["Read_length_median"] = round(np.median(lengths), 2)
        except Exception:
            summary["BWA_reads"] = "error"
            summary["Read_length_mean"] = "error"
            summary["Read_length_median"] = "error"

    # DamageProfiler summary
    dfile = os.path.join(args.damage_dir, f"damageprofiler_{pathogen_safe}", "5p_freq_misincorporations.txt")
    if os.path.exists(dfile):
        try:
            ddf = pd.read_csv(dfile, sep="\t", comment="#")
            if "C>T" in ddf.columns:
                summary["Damage_5p_CtoT"] = ddf["C>T"].iloc[0]
        except Exception:
            summary["Damage_5p_CtoT"] = "error"

    # Evenness and Coverage from Qualimap
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

    # Relative Entropy from file
    entropy_file = os.path.join(args.bam_dir, f"{args.sample}_{pathogen_safe}.mean_entropy.txt")
    if os.path.exists(entropy_file):
        try:
            with open(entropy_file) as f:
                summary["Relative_entropy"] = float(f.readline().strip())
        except Exception:
            summary["Relative_entropy"] = "error"

    # Expected Breadth from file
    expected_breadth_file = os.path.join(args.bam_dir, f"{args.sample}_{pathogen_safe}.expected_breadth.txt")
    if os.path.exists(expected_breadth_file):
        try:
            val = float(open(expected_breadth_file).readline().strip())
            summary["Expected_breadth"] = f"{round(val * 100, 2)}%"
        except Exception:
            summary["Expected_breadth"] = "error"

    # ANI from file
    ani_file = os.path.join(args.bam_dir, f"{args.sample}_{pathogen_safe}.ani.txt")
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

# -------------------- Save summary --------------------
os.makedirs(os.path.dirname(args.output), exist_ok=True)
df = pd.DataFrame(rows)

# Drop unused or deprecated columns
df.drop(columns=[
    "Entropy_plot",
    "Breadth_of_coverage",
    "Breadth_ratio",
    "Read_length"  # <- drop this old column if it still exists
], inplace=True, errors="ignore")

df.to_csv(args.output, index=False)
