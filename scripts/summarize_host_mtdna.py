import os
import pandas as pd
import pysam
import gzip
import sys
import re
import numpy as np

def extract_qualimap_coverage(report_path):
    """Extracts 'mean coverageData' from a Qualimap genome_results.txt file."""
    if not os.path.exists(report_path):
        return "NA"
    with open(report_path) as f:
        for line in f:
            match = re.search(r"mean coverageData\s*=\s*([\d\.]+)X", line)
            if match:
                return float(match.group(1))
    return "NA"

def read_contamination_warning(sample):
    """Read human contamination % and any other warnings from contamination log."""
    warning_file = f"logs/contamination_warnings/{sample}.txt"
    human_contam = "NA"
    other_warnings = []
    if os.path.exists(warning_file):
        with open(warning_file) as f:
            for line in f:
                line = line.strip()
                if line.lower().startswith("human contamination load"):
                    # Extract the float percentage from line
                    m = re.search(r"(\d+\.?\d*)%", line)
                    if m:
                        human_contam = float(m.group(1))
                elif "warning" in line.lower():
                    other_warnings.append(line)
    return human_contam, "; ".join(other_warnings) if other_warnings else "None"

def read_bam_read_lengths(bam_path):
    """Return mean and median read length from dedup bam."""
    if not os.path.exists(bam_path):
        return "NA", "NA"
    lengths = []
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    for read in bamfile.fetch(until_eof=True):
        if not read.is_unmapped:
            lengths.append(read.query_length)
    bamfile.close()
    if lengths:
        return round(np.mean(lengths), 2), round(np.median(lengths), 2)
    else:
        return "NA", "NA"

def read_damageprofiler_5pct(path):
    """Read mean 5pCtoT from DamageProfiler file."""
    if not os.path.exists(path):
        return "NA"
    values = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) == 2:
                try:
                    values.append(float(parts[1]))
                except ValueError:
                    continue
    if values:
        return round(np.mean(values), 4)
    else:
        return "NA"

if __name__ == "__main__":
    if "snakemake" in globals():
        samples = snakemake.params.samples.split(",")
        output_file = snakemake.output[0]
    else:
        samples = sys.argv[1].split(",")
        output_file = sys.argv[2]

    summary_data = []

    for sample in samples:
        # Species info
        species_file = f"results/{sample}/fastq_screen/{sample}_best_species.txt"
        species = "NA"
        if os.path.exists(species_file):
            with open(species_file) as f:
                species = f.read().strip()

        # Raw reads from config file
        r1_file = None
        with open("config/samples.tsv") as f:
            lines = f.readlines()[1:]  # skip header
            try:
                r1_file = [line.split("\t")[1] for line in lines if line.startswith(sample)][0]
            except IndexError:
                print(f"[WARNING] No R1 fastq entry found for sample '{sample}' in config.")

        raw_reads = 0
        if r1_file and os.path.exists(r1_file):
            with gzip.open(r1_file, "rt") as f:
                raw_reads = sum(1 for i, _ in enumerate(f) if i % 4 == 0)

        # Collapsed reads
        collapsed_file = f"results/{sample}/adapter_removal/{sample}.collapsed.gz"
        collapsed_reads = 0
        if os.path.exists(collapsed_file):
            with gzip.open(collapsed_file, "rt") as f:
                collapsed_reads = sum(1 for i, _ in enumerate(f) if i % 4 == 0)

        # Host Q30 BAM reads
        q30_bam_path = f"results/{sample}/bwa_host/{sample}_F4_q30.bam"
        q30_bai_path = q30_bam_path + ".bai"
        host_q30_reads = "NA"
        if os.path.exists(q30_bam_path) and os.path.exists(q30_bai_path):
            with pysam.AlignmentFile(q30_bam_path, "rb") as q30_bam:
                host_q30_reads = q30_bam.count()

        # Host deduplicated BAM reads
        dedup_bam_path = f"results/{sample}/bwa_host/{sample}.dedup.bam"
        dedup_bai_path = dedup_bam_path + ".bai"
        host_dedup_reads = "NA"
        if os.path.exists(dedup_bam_path) and os.path.exists(dedup_bai_path):
            with pysam.AlignmentFile(dedup_bam_path, "rb") as dedup_bam:
                host_dedup_reads = dedup_bam.count()

        # Host coverage
        host_qualimap_path = f"results/{sample}/qualimap/genome_results.txt"
        host_cov = extract_qualimap_coverage(host_qualimap_path)

        # mtDNA q30 BAM reads
        mt_q30_bam_path = f"results/{sample}/bwa_mtdna/{sample}_F4_q30.bam"
        mt_q30_bai_path = mt_q30_bam_path + ".bai"
        mt_q30_reads = "NA"
        if os.path.exists(mt_q30_bam_path) and os.path.exists(mt_q30_bai_path):
            with pysam.AlignmentFile(mt_q30_bam_path, "rb") as mt_q30_bam:
                mt_q30_reads = mt_q30_bam.count()

        # mtDNA deduplicated BAM reads
        mt_dedup_bam_path = f"results/{sample}/bwa_mtdna/{sample}.dedup.bam"
        mt_dedup_bai_path = mt_dedup_bam_path + ".bai"
        mt_dedup_reads = "NA"
        if os.path.exists(mt_dedup_bam_path) and os.path.exists(mt_dedup_bai_path):
            with pysam.AlignmentFile(mt_dedup_bam_path, "rb") as mt_dedup_bam:
                mt_dedup_reads = mt_dedup_bam.count()

        # mtDNA coverage
        mt_qualimap_path = f"results/{sample}/qualimap_mtdna/genome_results.txt"
        mt_cov = extract_qualimap_coverage(mt_qualimap_path)

        # Endogenous % for host
        try:
            host_endogenous = round((host_dedup_reads / raw_reads) * 100, 2) if isinstance(host_dedup_reads, int) and raw_reads > 0 else "NA"
        except:
            host_endogenous = "NA"

        # Duplication rate (host)
        duplication_rate = "NA"
        if isinstance(host_dedup_reads, int) and isinstance(host_q30_reads, int) and host_q30_reads > 0:
            duplication_rate = round(1 - (host_dedup_reads / host_q30_reads), 4)

        # Mean and median read length from dedup BAM (host)
        host_mean_len, host_median_len = read_bam_read_lengths(dedup_bam_path)

        # Mean and median read length from dedup BAM (mtDNA)
        mt_mean_len, mt_median_len = read_bam_read_lengths(mt_dedup_bam_path)

        # Human contamination and other warnings
        human_contam, other_warnings = read_contamination_warning(sample)

        # DamageProfiler 5pCtoT mean (host)
        damageprofiler_5pct_path = f"results/{sample}/damageprofiler_host/5pCtoT_freq.txt"
        damage_5pct_mean = read_damageprofiler_5pct(damageprofiler_5pct_path)

        summary_data.append({
            "sample": sample,
            "species": species,
            "raw_reads": raw_reads,
            "collapsed_reads": collapsed_reads,
            "host_q30_reads": host_q30_reads,
            "host_dedup_reads": host_dedup_reads,
            "host_coverage": host_cov,
            "host_endogenous_pct": host_endogenous,
            "duplication_rate": duplication_rate,
            "host_mean_read_length": host_mean_len,
            "host_median_read_length": host_median_len,
            "mt_q30_reads": mt_q30_reads,
            "mt_dedup_reads": mt_dedup_reads,
            "mt_coverage": mt_cov,
            "mt_mean_read_length": mt_mean_len,
            "mt_median_read_length": mt_median_len,
            "human_contamination_pct": human_contam,
            "damageprofiler_5pCtoT_mean": damage_5pct_mean,
            "other_warnings": other_warnings
        })

    df = pd.DataFrame(summary_data)
    df.to_excel(output_file, index=False)
