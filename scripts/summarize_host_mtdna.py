import os
import pandas as pd
import pysam
import gzip
import sys
import re
import numpy as np
import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

from sexing_utils import read_sexing_tsv, sexing_plot_path

def extract_qualimap_coverage(report_path):
    """Extracts 'mean coverageData' from a Qualimap genome_results.txt file."""
    if not os.path.exists(report_path):
        return "NA"
    try:
        with open(report_path) as f:
            for line in f:
                # Try multiple patterns to match Qualimap output format
                # Pattern 1: "mean coverageData = X.XX X"
                match = re.search(r"mean coverageData\s*=\s*([\d\.]+)\s*X", line, re.IGNORECASE)
                if match:
                    val = float(match.group(1))
                    return val if val > 0 else "NA"
                # Pattern 2: "mean coverage = X.XX"
                match = re.search(r"mean coverage\s*[=:]\s*([\d\.]+)", line, re.IGNORECASE)
                if match:
                    val = float(match.group(1))
                    return val if val > 0 else "NA"
    except Exception as e:
        print(f"[WARNING] Could not extract coverage from {report_path}: {e}")
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

def read_species_mismatch_warning(bio_sample):
    """Read species mismatch warning for a biological sample (sample-level, not PCR-level)."""
    host_warning_file = f"results/host/host_mapping/{bio_sample}.species_mismatch_warning.txt"
    mtdna_warning_file = f"results/host/mtdna_mapping/{bio_sample}.species_mismatch_warning.txt"
    warning_text = ""
    if os.path.exists(host_warning_file):
        with open(host_warning_file) as f:
            warning_text = f.read().strip()
    elif os.path.exists(mtdna_warning_file):
        with open(mtdna_warning_file) as f:
            warning_text = f.read().strip()
    return warning_text if warning_text else None

def read_bam_read_lengths(bam_path):
    """Return mean and median read length from dedup bam (sampled for speed)."""
    if not os.path.exists(bam_path):
        return "NA", "NA"
    lengths = []
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    # Sample up to 10000 reads for speed (or all if less than 10000)
    max_samples = 10000
    count = 0
    for read in bamfile.fetch(until_eof=True):
        if not read.is_unmapped:
            lengths.append(read.query_length)
            count += 1
            if count >= max_samples:
                break
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

def process_sample(sample, pcr_to_bio, samples_tsv_data):
    """Process a single sample and return summary data dictionary."""
    bio_sample = pcr_to_bio.get(sample, sample)
    # Species info
    species_file = f"results/libraries/{sample}/fastq_screen/{sample}_best_species.txt"
    species = "NA"
    if os.path.exists(species_file):
        with open(species_file) as f:
            species = f.read().strip()

    # Raw reads
    # Prefer pre-computed QC file (small text file from rule calculate_raw_reads).
    raw_reads_qc_file = f"results/libraries/{sample}/qc/{sample}.raw_reads.txt"
    raw_reads = 0

    if os.path.exists(raw_reads_qc_file):
        try:
            with open(raw_reads_qc_file) as f:
                line = f.readline().strip()
                if line:
                    raw_reads = int(float(line))
        except Exception as e:
            print(f"[WARNING] Could not read raw_reads from {raw_reads_qc_file}: {e}")
            raw_reads = 0
    else:
        # Fallback: compute directly from r1 FASTQ (legacy behaviour)
        r1_file = None
        for row in samples_tsv_data:
            pcr_id = row.get("pcr", row.get("sample", ""))
            if pcr_id == sample:
                r1_file = row.get("r1", row.get("R1", ""))
                break

        if r1_file and os.path.exists(r1_file):
            try:
                # Try using pigz if available (faster than zcat), otherwise use zcat
                # Count lines and divide by 4 (faster than grep for large files)
                decompressor = "pigz -dc" if subprocess.run(["which", "pigz"], capture_output=True).returncode == 0 else "zcat"
                result = subprocess.run(
                    f"{decompressor} {r1_file} | wc -l",
                    shell=True,
                    capture_output=True,
                    text=True,
                    timeout=600  # Increased timeout for very large files
                )
                if result.returncode == 0:
                    raw_reads = int(result.stdout.strip()) // 4
                else:
                    # Fallback: try with zcat if pigz failed
                    result = subprocess.run(
                        f"zcat {r1_file} | wc -l",
                        shell=True,
                        capture_output=True,
                        text=True,
                        timeout=600
                    )
                    if result.returncode == 0:
                        raw_reads = int(result.stdout.strip()) // 4
            except subprocess.TimeoutExpired:
                print(f"[WARNING] Timeout reading raw reads from {r1_file} (file may be very large)")
                raw_reads = 0
            except Exception as e:
                print(f"[WARNING] Could not read raw reads from {r1_file}: {e}")
                raw_reads = 0

    # Collapsed reads
    collapsed_reads = 0
    collapsed_reads_qc_file = f"results/libraries/{sample}/qc/{sample}.collapsed_reads.txt"
    collapsed_file = f"results/libraries/{sample}/adapter_removal/{sample}.collapsed.gz"

    # Prefer the precomputed collapsed_reads.txt (from calculate_collapsed_reads rule)
    if os.path.exists(collapsed_reads_qc_file):
        try:
            with open(collapsed_reads_qc_file) as f:
                line = f.readline().strip()
                if line:
                    collapsed_reads = int(float(line))
        except Exception as e:
            print(f"[WARNING] Could not read collapsed_reads from {collapsed_reads_qc_file}: {e}")
            collapsed_reads = 0
    # Backwards-compatible fallback: compute from the collapsed.gz if still present
    elif os.path.exists(collapsed_file):
        try:
            # Try using pigz if available (faster than zcat), otherwise use zcat
            # Count lines and divide by 4 (faster than grep for large files)
            decompressor = "pigz -dc" if subprocess.run(["which", "pigz"], capture_output=True).returncode == 0 else "zcat"
            result = subprocess.run(
                f"{decompressor} {collapsed_file} | wc -l",
                shell=True,
                capture_output=True,
                text=True,
                timeout=600  # Increased timeout for very large files
            )
            if result.returncode == 0:
                collapsed_reads = int(result.stdout.strip()) // 4
            else:
                # Fallback: try with zcat if pigz failed
                result = subprocess.run(
                    f"zcat {collapsed_file} | wc -l",
                    shell=True,
                    capture_output=True,
                    text=True,
                    timeout=600
                )
                if result.returncode == 0:
                    collapsed_reads = int(result.stdout.strip()) // 4
        except subprocess.TimeoutExpired:
            print(f"[WARNING] Timeout counting collapsed reads from {collapsed_file} (file may be very large)")
            collapsed_reads = 0
        except Exception as e:
            print(f"[WARNING] Could not count collapsed reads from {collapsed_file}: {e}")
            collapsed_reads = 0

    # Host deduplicated BAM reads (need this first for duplication rate calculation)
    dedup_bam_path = f"results/libraries/{sample}/host_mapping/{sample}.dedup.bam"
    dedup_bai_path = dedup_bam_path + ".bai"
    host_dedup_reads = "NA"
    if os.path.exists(dedup_bam_path):
        try:
            # Use samtools view -c for much faster counting
            result = subprocess.run(
                ["samtools", "view", "-c", "-F", "4", dedup_bam_path],
                capture_output=True,
                text=True,
                timeout=60
            )
            if result.returncode == 0:
                host_dedup_reads = int(result.stdout.strip())
        except Exception as e:
            print(f"[WARNING] Could not count reads from {dedup_bam_path}: {e}")
            host_dedup_reads = "NA"

    # Host Q30 reads and duplication rate from saved metrics file
    q30_metrics_file = f"results/libraries/{sample}/host_mapping/{sample}.q30_metrics.txt"
    host_q30_reads = "NA"
    duplication_rate = "NA"
    
    if os.path.exists(q30_metrics_file):
        try:
            with open(q30_metrics_file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("q30_reads="):
                        host_q30_reads = int(line.split("=")[1])
                    elif line.startswith("duplication_rate="):
                        dup_val = line.split("=")[1]
                        if dup_val != "NA":
                            duplication_rate = round(float(dup_val), 4)
        except Exception as e:
            print(f"[WARNING] Could not read q30 metrics from {q30_metrics_file}: {e}")
    
    # Fallback: try to read from BAM if metrics file doesn't exist
    if host_q30_reads == "NA":
        q30_bam_path = None
        if os.path.exists(f"results/libraries/{sample}/host_mapping/{sample}_F4_q30.sorted.bam"):
            q30_bam_path = f"results/libraries/{sample}/host_mapping/{sample}_F4_q30.sorted.bam"
        elif os.path.exists(f"results/libraries/{sample}/host_mapping/{sample}_F4_q30.bam"):
            q30_bam_path = f"results/libraries/{sample}/host_mapping/{sample}_F4_q30.bam"
        
        if q30_bam_path and os.path.exists(q30_bam_path):
            try:
                # Use samtools view -c for much faster counting
                result = subprocess.run(
                    ["samtools", "view", "-c", "-F", "4", q30_bam_path],
                    capture_output=True,
                    text=True,
                    timeout=60
                )
                if result.returncode == 0:
                    host_q30_reads = int(result.stdout.strip())
                    # Calculate duplication rate if we have both counts
                    if isinstance(host_dedup_reads, int) and isinstance(host_q30_reads, int) and host_q30_reads > 0:
                        duplication_rate = round(1 - (host_dedup_reads / host_q30_reads), 4)
            except:
                host_q30_reads = "NA"

    # Host coverage
    host_qualimap_path = f"results/libraries/{sample}/qualimap/genome_results.txt"
    host_cov = extract_qualimap_coverage(host_qualimap_path)

    # mtDNA Q30 reads from saved metrics file
    mt_q30_metrics_file = f"results/libraries/{sample}/mtdna_mapping/{sample}.q30_metrics.txt"
    mt_q30_reads = "NA"
    
    if os.path.exists(mt_q30_metrics_file):
        try:
            with open(mt_q30_metrics_file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("q30_reads="):
                        mt_q30_reads = int(line.split("=")[1])
        except Exception as e:
            print(f"[WARNING] Could not read mtDNA q30 metrics from {mt_q30_metrics_file}: {e}")
    
    # Fallback: try to read from BAM if metrics file doesn't exist
    if mt_q30_reads == "NA":
        mt_q30_bam_path = None
        if os.path.exists(f"results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.sorted.bam"):
            mt_q30_bam_path = f"results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.sorted.bam"
        elif os.path.exists(f"results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.bam"):
            mt_q30_bam_path = f"results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.bam"
        
        if mt_q30_bam_path and os.path.exists(mt_q30_bam_path):
            try:
                # Use samtools view -c for much faster counting
                result = subprocess.run(
                    ["samtools", "view", "-c", "-F", "4", mt_q30_bam_path],
                    capture_output=True,
                    text=True,
                    timeout=60
                )
                if result.returncode == 0:
                    mt_q30_reads = int(result.stdout.strip())
            except:
                mt_q30_reads = "NA"

    # mtDNA deduplicated BAM reads
    mt_dedup_bam_path = f"results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam"
    mt_dedup_bai_path = mt_dedup_bam_path + ".bai"
    mt_dedup_reads = "NA"
    if os.path.exists(mt_dedup_bam_path):
        try:
            # Use samtools view -c for much faster counting
            result = subprocess.run(
                ["samtools", "view", "-c", "-F", "4", mt_dedup_bam_path],
                capture_output=True,
                text=True,
                timeout=60
            )
            if result.returncode == 0:
                mt_dedup_reads = int(result.stdout.strip())
        except Exception as e:
            print(f"[WARNING] Could not count reads from {mt_dedup_bam_path}: {e}")
            mt_dedup_reads = "NA"

    # mtDNA coverage
    mt_qualimap_path = f"results/libraries/{sample}/qualimap_mtdna/genome_results.txt"
    mt_cov = extract_qualimap_coverage(mt_qualimap_path)

    # Endogenous % for host (PCR level: dedup reads / raw reads * 100)
    try:
        host_endogenous = round((host_dedup_reads / raw_reads) * 100, 2) if isinstance(host_dedup_reads, int) and raw_reads > 0 else "NA"
    except:
        host_endogenous = "NA"
    
    # Duplication rate already read from metrics file above, if not available calculate it
    if duplication_rate == "NA" and isinstance(host_dedup_reads, int) and isinstance(host_q30_reads, int) and host_q30_reads > 0:
        duplication_rate = round(1 - (host_dedup_reads / host_q30_reads), 4)

    # Mean and median read length from dedup BAM (host)
    host_mean_len, host_median_len = read_bam_read_lengths(dedup_bam_path)

    # Mean and median read length from dedup BAM (mtDNA)
    mt_mean_len, mt_median_len = read_bam_read_lengths(mt_dedup_bam_path)

    # Human contamination and other warnings
    human_contam, other_warnings = read_contamination_warning(sample)

    # DamageProfiler 5pCtoT mean (host)
    damageprofiler_5pct_path = f"results/libraries/{sample}/damageprofiler_host/5pCtoT_freq.txt"
    damage_5pct_mean = read_damageprofiler_5pct(damageprofiler_5pct_path)

    return {
        "pcr_id": sample,
        "bio_sample": bio_sample,
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
        "mtdna_coverage": mt_cov,  # Changed from mt_coverage to mtdna_coverage for consistency
        "mt_mean_read_length": mt_mean_len,
        "mt_median_read_length": mt_median_len,
        "human_contamination_pct": human_contam,
        "damageprofiler_5pCtoT_mean": damage_5pct_mean,
        "other_warnings": other_warnings
    }

if __name__ == "__main__":
    if "snakemake" in globals():
        samples = snakemake.params.samples.split(",")  # these are PCR IDs
        output_file = snakemake.output[0]
        samples_tsv = snakemake.input.samples_tsv
        cleanup_enabled = snakemake.config.get("cleanup_intermediates", False)
        # Use threads allocated by Snakemake
        available_threads = snakemake.threads
    else:
        samples = sys.argv[1].split(",")
        output_file = sys.argv[2]
        samples_tsv = "config/samples.tsv"
        cleanup_enabled = False
        # Fallback to CPU count when not running under Snakemake
        available_threads = os.cpu_count() or 8

    # Map PCR IDs to biological samples using the samples.tsv
    pcr_to_bio = {}
    # Pre-load samples_tsv data to avoid reading it multiple times
    samples_tsv_data = []
    with open(samples_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        # Normalize header names to avoid issues with BOMs or stray whitespace
        if reader.fieldnames is not None:
            reader.fieldnames = [fn.strip().lstrip("\ufeff") for fn in reader.fieldnames]
        has_pcr_col = "pcr" in reader.fieldnames if reader.fieldnames is not None else False
        for row in reader:
            # Skip completely empty / malformed rows
            if not row or "sample" not in row or row["sample"] is None or not str(row["sample"]).strip():
                continue
            bio = str(row["sample"]).strip()
            pcr = str(row["pcr"]).strip() if has_pcr_col and row.get("pcr") else bio
            pcr_to_bio[pcr] = bio
            samples_tsv_data.append(row)

    # Process samples in parallel using ThreadPoolExecutor
    # Use threads allocated by Snakemake (or CPU count if not running under Snakemake)
    # Cap at number of samples to avoid creating unnecessary threads
    max_workers = min(available_threads, len(samples))
    print(f"[INFO] Processing {len(samples)} samples with {max_workers} threads (available: {available_threads})...")
    
    summary_data = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_sample = {
            executor.submit(process_sample, sample, pcr_to_bio, samples_tsv_data): sample
            for sample in samples
        }
        
        # Collect results as they complete
        completed = 0
        for future in as_completed(future_to_sample):
            sample = future_to_sample[future]
            try:
                result = future.result()
                summary_data.append(result)
                completed += 1
                if completed % 10 == 0:
                    print(f"[INFO] Processed {completed}/{len(samples)} samples...")
            except Exception as e:
                print(f"[ERROR] Failed to process sample {sample}: {e}")
                # Add a placeholder entry with NA values
                summary_data.append({
                    "pcr_id": sample,
                    "bio_sample": pcr_to_bio.get(sample, sample),
                    "species": "NA",
                    "raw_reads": "NA",
                    "collapsed_reads": "NA",
                    "host_q30_reads": "NA",
                    "host_dedup_reads": "NA",
                    "host_coverage": "NA",
                    "host_endogenous_pct": "NA",
                    "duplication_rate": "NA",
                    "host_mean_read_length": "NA",
                    "host_median_read_length": "NA",
                    "mt_q30_reads": "NA",
                    "mt_dedup_reads": "NA",
                    "mtdna_coverage": "NA",  # Changed from mt_coverage to mtdna_coverage for consistency
                    "mt_mean_read_length": "NA",
                    "mt_median_read_length": "NA",
                    "human_contamination_pct": "NA",
                    "damageprofiler_5pCtoT_mean": "NA",
                    "other_warnings": "NA",
                })
    
    print(f"[INFO] Completed processing all {len(samples)} samples.")

    # PCR-level DataFrame
    pcr_df = pd.DataFrame(summary_data)

    # Build sample-level aggregates
    sample_rows = []
    for bio_sample, grp in pcr_df.groupby("bio_sample"):
        row = {"bio_sample": bio_sample}

        # Take first non-NA species if available
        species_vals = [v for v in grp["species"].tolist() if v != "NA"]
        row["species"] = species_vals[0] if species_vals else "NA"

        sexing_tsv = f"results/samples/{bio_sample}/sexing/{bio_sample}_sexing.tsv"
        sexing = read_sexing_tsv(sexing_tsv)
        sexing["sexing_plot"] = sexing_plot_path(bio_sample)
        row.update(sexing)

        def sum_numeric(col):
            vals = pd.to_numeric(grp[col], errors="coerce")
            s = vals.sum(skipna=True)
            return s if not pd.isna(s) else "NA"

        def mean_numeric(col):
            vals = pd.to_numeric(grp[col], errors="coerce")
            m = vals.mean(skipna=True)
            return round(m, 4) if not pd.isna(m) else "NA"

        # Summed counts
        for col in [
            "raw_reads",
            "collapsed_reads",
            "host_q30_reads",
            "host_dedup_reads",
            "mt_q30_reads",
            "mt_dedup_reads",
        ]:
            if col in grp.columns:
                row[col] = sum_numeric(col)

        # Sample-level coverage from merged BAM Qualimap (not averaged from PCRs)
        host_sample_qm = f"results/samples/{bio_sample}/qualimap/genome_results.txt"
        mt_sample_qm = f"results/samples/{bio_sample}/qualimap_mtdna/genome_results.txt"
        row["host_coverage"] = extract_qualimap_coverage(host_sample_qm)
        row["mtdna_coverage"] = extract_qualimap_coverage(mt_sample_qm)  # Changed from mt_coverage to mtdna_coverage for consistency

        # Sample-level endogenous % for host: (sum of host_dedup_reads / sum of raw_reads) * 100
        try:
            host_dedup_sum = row.get("host_dedup_reads")
            raw_reads_sum = row.get("raw_reads")
            # Convert to numeric, handling both numeric types and "NA" strings
            if host_dedup_sum != "NA" and raw_reads_sum != "NA":
                host_dedup_sum = pd.to_numeric(host_dedup_sum, errors="coerce")
                raw_reads_sum = pd.to_numeric(raw_reads_sum, errors="coerce")
                if not pd.isna(host_dedup_sum) and not pd.isna(raw_reads_sum) and raw_reads_sum > 0:
                    row["host_endogenous_pct"] = round((host_dedup_sum / raw_reads_sum) * 100, 2)
                else:
                    row["host_endogenous_pct"] = "NA"
            else:
                row["host_endogenous_pct"] = "NA"
        except Exception as e:
            print(f"[WARNING] Could not calculate sample-level endogenous for {bio_sample}: {e}")
            row["host_endogenous_pct"] = "NA"
        
        # Sample-level duplication rate: (1 - host_dedup_reads / host_q30_reads) * 100
        # Note: host_q30_reads and host_dedup_reads are already set in the loop above
        try:
            host_q30_sum = row.get("host_q30_reads")
            host_dedup_sum = row.get("host_dedup_reads")
            # Convert to numeric, handling both numeric types and "NA" strings
            if host_q30_sum != "NA" and host_dedup_sum != "NA":
                host_q30_sum = pd.to_numeric(host_q30_sum, errors="coerce")
                host_dedup_sum = pd.to_numeric(host_dedup_sum, errors="coerce")
                if not pd.isna(host_q30_sum) and not pd.isna(host_dedup_sum) and host_q30_sum > 0:
                    row["duplication_rate"] = round((1 - (host_dedup_sum / host_q30_sum)) * 100, 2)
                else:
                    row["duplication_rate"] = "NA"
            else:
                row["duplication_rate"] = "NA"
        except Exception as e:
            print(f"[WARNING] Could not calculate sample-level duplication rate for {bio_sample}: {e}")
            row["duplication_rate"] = "NA"

        # Mean read lengths & damage (still averaged over PCRs)
        for col in [
            "host_mean_read_length",
            "host_median_read_length",
            "mt_mean_read_length",
            "mt_median_read_length",
            "damageprofiler_5pCtoT_mean",
        ]:
            if col in grp.columns:
                row[col] = mean_numeric(col)

        # Contamination: take max human_contamination_pct, and concatenate warnings
        if "human_contamination_pct" in grp.columns:
            vals = pd.to_numeric(grp["human_contamination_pct"], errors="coerce")
            mx = vals.max(skipna=True)
            row["human_contamination_pct"] = round(mx, 4) if not pd.isna(mx) else "NA"

        if "other_warnings" in grp.columns:
            warnings = set(
                w for w in grp["other_warnings"].tolist()
                if isinstance(w, str) and w != "None"
            )
            row["other_warnings"] = "; ".join(sorted(warnings)) if warnings else "None"

        # Check for species mismatch warning (sample-level)
        species_mismatch_warning = read_species_mismatch_warning(bio_sample)
        if species_mismatch_warning:
            # Add species mismatch warning to other_warnings
            if row.get("other_warnings", "None") == "None":
                row["other_warnings"] = f"SPECIES_MISMATCH: {species_mismatch_warning}"
            else:
                row["other_warnings"] = f"{row['other_warnings']}; SPECIES_MISMATCH: {species_mismatch_warning}"
            # Set host/mtdna metrics to "NA" or warning message when species mismatch detected
            row["host_coverage"] = "SPECIES_MISMATCH"
            row["mtdna_coverage"] = "SPECIES_MISMATCH"
            row["host_dedup_reads"] = "SPECIES_MISMATCH"
            row["mt_dedup_reads"] = "SPECIES_MISMATCH"
            row["host_q30_reads"] = "SPECIES_MISMATCH"
            row["mt_q30_reads"] = "SPECIES_MISMATCH"
            row["host_endogenous_pct"] = "SPECIES_MISMATCH"
            row["duplication_rate"] = "SPECIES_MISMATCH"
            # Keep raw_reads and collapsed_reads as they are still valid

        sample_rows.append(row)

    sample_df = pd.DataFrame(sample_rows)

    # Write both levels to a single Excel file
    sexing_cols = [
        "bio_sample",
        "species",
        "sexing_call",
        "sexing_female_prob",
        "sexing_likelihood_ratio",
        "sexing_status",
        "sexing_plot",
        "sexing_note",
    ]
    sexing_sample_df = sample_df[[c for c in sexing_cols if c in sample_df.columns]].copy()

    with pd.ExcelWriter(output_file) as writer:
        pcr_df.to_excel(writer, sheet_name="PCR_level", index=False)
        sample_df.to_excel(writer, sheet_name="Sample_level", index=False)
        if not sexing_sample_df.empty:
            sexing_sample_df.to_excel(writer, sheet_name="Sexing", index=False)

    # Cleanup intermediate files if enabled
    if cleanup_enabled:
        import glob
        import shutil
        print("[CLEANUP] Removing intermediate files...")
        for sample in samples:
            # Remove intermediate host BAMs
            intermediate_host_files = [
                f"results/libraries/{sample}/host_mapping/{sample}_F4.bam",
                f"results/libraries/{sample}/host_mapping/{sample}_F4.bam.bai",
                f"results/libraries/{sample}/host_mapping/{sample}_F4_q30.bam",
                f"results/libraries/{sample}/host_mapping/{sample}_F4_q30.bam.bai",
                f"results/libraries/{sample}/host_mapping/{sample}_F4_q30.sorted.bam",
                f"results/libraries/{sample}/host_mapping/{sample}_F4.sai",
            ]
            # Remove intermediate mtDNA BAMs
            intermediate_mtdna_files = [
                f"results/libraries/{sample}/mtdna_mapping/{sample}_F4.bam",
                f"results/libraries/{sample}/mtdna_mapping/{sample}_F4.bam.bai",
                f"results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.bam",
                f"results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.bam.bai",
                f"results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.sorted.bam",
                f"results/libraries/{sample}/mtdna_mapping/{sample}_F4.sai",
            ]
            # Remove prinseq uncompressed files
            prinseq_files = [
                f"results/libraries/{sample}/prinseq/{sample}-passed.fq",
            ]
            # Remove empty host-unmapped files (if they're just placeholders)
            host_unmapped_files = [
                f"results/libraries/{sample}/host_mapping/{sample}_host_unaligned.bam",
                f"results/libraries/{sample}/host_mapping/{sample}_host_unaligned.fastq.gz",
            ]
            
            all_intermediate = (intermediate_host_files + intermediate_mtdna_files + 
                             prinseq_files + host_unmapped_files)
            
            for filepath in all_intermediate:
                if os.path.exists(filepath):
                    # For empty FASTQ files (placeholders), check size first
                    if filepath.endswith(".fastq.gz"):
                        try:
                            if os.path.getsize(filepath) == 0:
                                os.remove(filepath)
                                print(f"  Removed empty placeholder: {filepath}")
                        except:
                            pass
                    else:
                        try:
                            os.remove(filepath)
                            print(f"  Removed: {filepath}")
                        except Exception as e:
                            print(f"  Warning: Could not remove {filepath}: {e}")
        
        print("[CLEANUP] Intermediate file cleanup complete.")
