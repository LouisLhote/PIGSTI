#!/usr/bin/env python3
"""
Diagnostic script to check how many samples PIGSTI will process.
Run this from your PIGSTI directory to see:
- Which samples.tsv file is being used
- How many samples/PCRs are loaded
- Any issues with sample loading
"""

import sys
import os
import csv
import yaml

# Add PIGSTI directory to path if needed
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Load config
config_file = "config/config.yaml"
if not os.path.exists(config_file):
    print(f"ERROR: Config file not found: {config_file}")
    print(f"Current directory: {os.getcwd()}")
    sys.exit(1)

with open(config_file) as f:
    CFG = yaml.safe_load(f)

SAMPLES_TSV = CFG.get("samples", "config/samples.tsv")
print(f"Config file: {config_file}")
print(f"Samples TSV path: {SAMPLES_TSV}")
print(f"Full path: {os.path.abspath(SAMPLES_TSV)}")
print(f"File exists: {os.path.exists(SAMPLES_TSV)}")
print()

if not os.path.exists(SAMPLES_TSV):
    print(f"ERROR: Samples TSV not found: {SAMPLES_TSV}")
    sys.exit(1)

# Load samples like Snakefile does
SAMPLES_DICT = {}
PCR_INFO = {}
SAMPLE_TO_PCRS = {}
PCRS = []

with open(SAMPLES_TSV) as f:
    reader = csv.DictReader(f, delimiter="\t")
    # Normalize header names to avoid issues with BOMs or stray whitespace
    if reader.fieldnames is not None:
        reader.fieldnames = [fn.strip().lstrip("\ufeff") for fn in reader.fieldnames]
    has_pcr_col = "pcr" in reader.fieldnames if reader.fieldnames is not None else False
    
    print(f"Columns in TSV: {reader.fieldnames}")
    print(f"Has 'pcr' column: {has_pcr_col}")
    print()
    
    row_count = 0
    skipped_count = 0
    
    for row in reader:
        row_count += 1
        # Skip completely empty / malformed rows
        if not row or "sample" not in row or row["sample"] is None or not str(row["sample"]).strip():
            skipped_count += 1
            print(f"  [SKIP] Row {row_count}: Empty or missing 'sample' column")
            continue

        sample = str(row["sample"]).strip()
        pcr_id = row["pcr"] if has_pcr_col and row.get("pcr") else sample
        pcr_id = str(pcr_id).strip() if pcr_id else sample

        # Per-PCR read lists
        if pcr_id not in SAMPLES_DICT:
            SAMPLES_DICT[pcr_id] = {"sample": sample, "r1": [], "r2": []}
        SAMPLES_DICT[pcr_id]["r1"].append(row["r1"])
        if row["r2"] and row["r2"].strip():
            SAMPLES_DICT[pcr_id]["r2"].append(row["r2"])

        # PCR metadata
        if pcr_id not in PCR_INFO:
            source = row.get("source", "").strip() if row.get("source") is not None else ""
            if not source:
                source = "LOCAL"
            PCR_INFO[pcr_id] = {
                "sample": sample,
                "r1": row["r1"],
                "r2": row["r2"],
                "RGLB": row.get("RGLB", ""),
                "sequencing_run": row.get("sequencing_run", ""),
                "source": source,
            }
            PCRS.append(pcr_id)
            if sample not in SAMPLE_TO_PCRS:
                SAMPLE_TO_PCRS[sample] = []
            SAMPLE_TO_PCRS[sample].append(pcr_id)

SAMPLES = list(SAMPLES_DICT.keys())
BIO_SAMPLES = list(SAMPLE_TO_PCRS.keys())

print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Total rows in TSV: {row_count}")
print(f"Skipped rows: {skipped_count}")
print(f"Loaded PCRs: {len(SAMPLES)}")
print(f"Biological samples: {len(BIO_SAMPLES)}")
print()
print("First 10 PCRs:")
for pcr in SAMPLES[:10]:
    print(f"  {pcr} -> {PCR_INFO[pcr]['sample']}")
if len(SAMPLES) > 10:
    print(f"  ... and {len(SAMPLES) - 10} more")
print()
print("First 10 biological samples:")
for bio in BIO_SAMPLES[:10]:
    pcr_list = SAMPLE_TO_PCRS[bio]
    print(f"  {bio} -> {len(pcr_list)} PCR(s): {', '.join(pcr_list[:3])}{'...' if len(pcr_list) > 3 else ''}")
if len(BIO_SAMPLES) > 10:
    print(f"  ... and {len(BIO_SAMPLES) - 10} more")
print()
print("=" * 60)
print("This is what Snakemake will process:")
print(f"  - {len(SAMPLES)} PCR-level rules (adapter removal, host/mtDNA mapping per PCR)")
print(f"  - {len(BIO_SAMPLES)} sample-level rules (KrakenUniq, HOPS, summaries, merged BAMs)")
print("=" * 60)


