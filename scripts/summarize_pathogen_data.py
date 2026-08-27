#!/usr/bin/env python

import os
import sys
from pathlib import Path

import pandas as pd
import argparse
import pysam
import numpy as np
import glob
import json
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pigsti_naming import safe_pathogen_name, hops_species_token

# -------------------- Argument parsing --------------------
parser = argparse.ArgumentParser()
parser.add_argument("--sample", required=True)
parser.add_argument(
    "--evalue",
    "--escore",
    dest="evalue",
    required=True,
    help="E-value pathogen CSV (results/pathogen/{bio}/evalue/pathogen/...)",
)
parser.add_argument("--hops", required=True)
parser.add_argument("--spreadsheet", required=True)
parser.add_argument("--bam_dir", required=True)
parser.add_argument("--qualimap_dir", required=True)
parser.add_argument("--damage_dir", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--cleanup", action="store_true", help="Clean up intermediate files after collecting data")
parser.add_argument("--config", default="config/config.yaml", help="Main pipeline config (for pathogen_detection_criteria)")
args = parser.parse_args()

# -------------------- Load detection criteria config --------------------
_default_criteria = {
    "ani_threshold": 96.5,
    "entropy_threshold": 0.9,
    "entropy_threshold_virus": 0.7,
    "breadth_ratio_threshold": 0.8,
    "damage_5p_ct_threshold": 0.01,
    "damage_3p_ga_threshold": 0.01,
    "hops_detection_threshold": 2,
    "hops_edit_distance_threshold": 3,
    "hops_damage_threshold": 4,
    "high_confidence_threshold": 0.7,
    "edit_distance_decay_quality_threshold": 0.65,
    "edit_distance_logr2_threshold": 0.60,
    "read_mapping_ratio_threshold": 0.5,
    "escore_threshold": 5,
    "reads_threshold": 50,
    "guellil_evalue_threshold": 0.001,
    "use_evalue_for_detection": True,
    # edit_distance_no_damage: prefer *_min (>=) in config.yaml; legacy *_max (<=) supported.
    "edit_distance_damage_minus_no_damage_min": 0.10,
}

criteria_config = dict(_default_criteria)

# Prefer criteria embedded in main config.yaml (also reuse this config for metadata fields)
pipeline_cfg = {}
if args.config and os.path.exists(args.config):
    try:
        with open(args.config, "r") as f:
            pipeline_cfg = yaml.safe_load(f) or {}
    except Exception:
        pipeline_cfg = {}

if isinstance(pipeline_cfg, dict) and isinstance(pipeline_cfg.get("pathogen_detection_criteria"), dict):
    criteria_config.update(pipeline_cfg["pathogen_detection_criteria"])

# Backwards compatible fallback: YAML next to spreadsheet
criteria_config_path = os.path.join(os.path.dirname(args.spreadsheet), "pathogen_detection_criteria.yaml")
if os.path.exists(criteria_config_path):
    try:
        with open(criteria_config_path, "r") as f:
            legacy = yaml.safe_load(f) or {}
        if isinstance(legacy, dict) and not (args.config and os.path.exists(args.config)):
            # If config.yaml exists we already loaded it; still allow legacy to fill missing keys.
            criteria_config.update(legacy)
        elif isinstance(legacy, dict):
            # fill missing keys only
            for k, v in legacy.items():
                if k not in criteria_config:
                    criteria_config[k] = v
    except Exception:
        pass

# -------------------- Host / sample metadata (best-guess, screening-only compatible) --------------------
# In screening-only mode there may be no host BAM/QC, but fastq_screen still produces best_species.txt.
host_best_species = "NA"
_candidate_species_files = []

# 1) Sample-level best species (exists in some pipeline modes)
_candidate_species_files.append(
    os.path.join("results", args.sample, "fastq_screen", f"{args.sample}_best_species.txt")
)

# 2) PCR-level best species: in many setups fastq_screen is run per-PCR, not per-bio-sample.
#    Use config/samples.tsv to map bio-sample -> PCR IDs, then read the first available.
try:
    samples_tsv = None
    if isinstance(pipeline_cfg, dict) and pipeline_cfg.get("samples"):
        samples_tsv = str(pipeline_cfg.get("samples"))
    if samples_tsv and os.path.exists(samples_tsv):
        st = pd.read_csv(samples_tsv, sep="\t")
        if "sample" in st.columns and "pcr" in st.columns:
            pcrs = st.loc[st["sample"].astype(str) == str(args.sample), "pcr"].dropna().astype(str).tolist()
            for pcr in pcrs:
                _candidate_species_files.append(
                    os.path.join("results", pcr, "fastq_screen", f"{pcr}_best_species.txt")
                )
except Exception:
    # Best-effort only; keep host_best_species = "NA" if we can't map.
    pass

for best_species_path in _candidate_species_files:
    if os.path.exists(best_species_path):
        try:
            with open(best_species_path, "r", encoding="utf-8", errors="ignore") as f:
                host_best_species = (f.readline() or "").strip() or "NA"
            if host_best_species != "NA":
                break
        except Exception:
            host_best_species = "NA"

host_aligner = "NA"
pathogen_screening_only = False
if isinstance(pipeline_cfg, dict):
    host_aligner = str(pipeline_cfg.get("host_aligner", "NA") or "NA")
    pathogen_screening_only = bool(pipeline_cfg.get("pathogen_screening_only", False))

# -------------------- Load reference spreadsheet --------------------
spreadsheet = pd.read_csv(args.spreadsheet)
spreadsheet["Krakenuniq name"] = spreadsheet["Krakenuniq name"].str.strip()
spreadsheet["Hops name"] = spreadsheet["Hops name"].str.strip()

# -------------------- Load escore (Krakenuniq) --------------------
escore_df = pd.read_csv(args.evalue)
escore_df.columns = escore_df.columns.str.strip()
escore_df["taxonomy"] = escore_df["taxonomy"].str.strip()

# Decide which KrakenUniq read-count column to use
# KrakenUniq report header is typically:
# "%\treads\ttaxReads\tkmers\tdup\tcov\ttaxID\trank\ttaxName"
#  - 'reads'    = reads in clade (including children)
#  - 'taxReads' = reads assigned directly to this taxon
# For detection thresholds we want the clade-level read count -> 'reads' when available.
if "reads" in escore_df.columns:
    reads_col = "reads"
elif "taxReads" in escore_df.columns:
    reads_col = "taxReads"
elif "# of reads" in escore_df.columns:
    reads_col = "# of reads"
else:
    raise ValueError(
        "Missing read count column in Escore file (expected 'reads', 'taxReads', or '# of reads')"
    )

# -------------------- Load HOPS (optional) --------------------
if args.hops and os.path.exists(args.hops):
    hops_df = pd.read_csv(args.hops, sep="\t")
    hops_df.columns = hops_df.columns.str.replace('"', '').str.strip()
    hops_df.rename(columns={"node": "Species"}, inplace=True)
    hops_df["Species"] = hops_df["Species"].str.replace('"', '').str.strip()
    sample_col = f"{args.sample}_unaligned.rma6"
else:
    # HOPS disabled or file missing: use an empty dataframe and skip HOPS-based metrics
    hops_df = pd.DataFrame()
    sample_col = None

# -------------------- Summarize detected pathogens only --------------------
rows = []

detected_kraken = set(escore_df["taxonomy"])
detected_hops = set()
if sample_col is not None and not hops_df.empty and sample_col in hops_df.columns:
    detected_hops = set(hops_df[hops_df[sample_col] > 1]["Species"])

# Apply the same filtering logic as get_sample_ref_pairs() in Snakefile
# Only include pathogens that pass E-value/Escore + reads thresholds
observed = set()
use_evalue = bool(criteria_config.get("use_evalue_for_detection", False))
escore_min_default = criteria_config.get("escore_threshold", 5)
reads_min_default = criteria_config.get("reads_threshold", 50)
evalue_max_default = criteria_config.get("guellil_evalue_threshold", 0.001)

for _, row in spreadsheet.iterrows():
    kraken_name = row["Krakenuniq name"]
    
    # Check if detected by KrakenUniq
    if kraken_name not in detected_kraken:
        # Also check HOPS
        if hops_species_token(row["Hops name"]) in detected_hops:
            observed.add(kraken_name)
        continue
    
    # Get pathogen-specific thresholds
    pathogen_row = spreadsheet[spreadsheet["Krakenuniq name"] == kraken_name]
    min_escore = escore_min_default
    min_reads = reads_min_default
    max_evalue = evalue_max_default
    
    if not pathogen_row.empty:
        if "min_escore" in pathogen_row.columns:
            val = pathogen_row.iloc[0]["min_escore"]
            if pd.notna(val):
                try:
                    min_escore = float(val)
                except Exception:
                    pass
        if "min_reads" in pathogen_row.columns:
            val = pathogen_row.iloc[0]["min_reads"]
            if pd.notna(val):
                try:
                    min_reads = int(val)
                except Exception:
                    pass
        # E-value threshold
        for colname in ["Guellil_et_al_Evalue_threshold", "max_evalue", "evalue_threshold"]:
            if colname in pathogen_row.columns:
                val = pathogen_row.iloc[0][colname]
                if pd.notna(val):
                    try:
                        max_evalue = float(val)
                    except Exception:
                        pass
                break
    
    # Get Escore data for this pathogen
    match = escore_df[escore_df["taxonomy"] == kraken_name]
    if match.empty:
        continue
    
    # Get reads count
    reads_val = 0
    if reads_col in match.columns:
        try:
            reads_val = int(match[reads_col].values[0])
        except Exception:
            reads_val = 0
    
    reads_ok = reads_val >= min_reads
    if not reads_ok:
        continue
    
    # Check E-value or Escore threshold
    if use_evalue:
        # E-value mode: require E-value column and check threshold
        if "Guellil_et_al_Evalue" in escore_df.columns:
            try:
                evalue_raw = match["Guellil_et_al_Evalue"].values[0]
                if pd.notna(evalue_raw):
                    # Convert to float (handles scientific notation like "4.51E-05")
                    evalue_val = float(evalue_raw)
                    # E-value metric must be > threshold (larger is better)
                    if evalue_val > max_evalue:
                        observed.add(kraken_name)
                    # If evalue_val <= max_evalue, do NOT add (explicitly skip)
            except (ValueError, TypeError, KeyError, IndexError):
                pass
        # If E-value mode enabled but column missing, skip this pathogen (don't fall back to Escore)
    else:
        # Escore mode (larger is better, must be >= threshold)
        try:
            escore_raw = match["Escore"].values[0]
            if pd.notna(escore_raw):
                escore_val = float(escore_raw)
                if escore_val >= min_escore:
                    observed.add(kraken_name)
        except (ValueError, TypeError, KeyError, IndexError):
            pass

# If nothing passed strict filtering, still emit a useful per-sample summary.
# This is especially important in screening-only mode where users still want
# the basic read counts/E-scores/E-values in the spreadsheet even if all
# detections fail thresholds.
include_threshold_failures = bool(
    (isinstance(pipeline_cfg, dict) and pipeline_cfg.get("pathogen_summary_include_threshold_failures", False))
    or pathogen_screening_only
)
if include_threshold_failures and len(observed) == 0:
    for _, row in spreadsheet.iterrows():
        kraken_name = row["Krakenuniq name"]
        if kraken_name in detected_kraken:
            observed.add(kraken_name)

for _, pathogen_row in spreadsheet.iterrows():
    kraken_name = pathogen_row["Krakenuniq name"]
    if kraken_name not in observed:
        continue

    hops_name = pathogen_row["Hops name"]
    pathogen_safe = safe_pathogen_name(kraken_name)
    hops_name_safe = hops_species_token(hops_name)

    summary = {
        "Sample": args.sample,
        "Pathogen": kraken_name,
        "Host_best_species": host_best_species,
        "Host_aligner": host_aligner,
        "Pathogen_screening_only": pathogen_screening_only,
        "Krakenuniq_reads": 0,
        "Escore": "NA",
        "Guellil_et_al_Evalue": "NA",
        "Detected_by_Krakenuniq": False,
        "HOPS_score": 0,
        "Detected_by_HOPS": False,
        "Mapped_reads": "NA",  # Changed from BWA_reads to Mapped_reads (works for both BWA and Bowtie2)
        "Read_mapping_ratio": "NA",  # Mapped_reads / Krakenuniq_reads (consistency metric)
        "Coverage": "NA",
        "Evenness": "NA",
        "Read_length_mean": "NA",
        "Read_length_median": "NA",
        "Damage_5p_CtoT": "NA",
        "Relative_entropy": "NA",
        "Breadth_ratio": "NA",  # Ratio of observed to expected breadth (as in Nature paper)
        "ANI": "NA",
        "Edit_distance_decay_quality": "NA",
        "Edit_distance_decay_quality_default": "NA",
        "Genus_ranking": "NA",  # Ranking position within genus
        "Score": "NA",  # Combined score format: "X/Y" (e.g., "3/6" or "5/9")
        "Criteria_passed": "NA"  # Comma-separated list of passed criteria
        # Removed Entropy_plot, Breadth_of_coverage, Expected_breadth columns
        # Removed Detection_score, Max_possible_score, Criteria_breakdown (replaced with Score and Criteria_passed)
    }

    # Krakenuniq match
    match = escore_df[escore_df["taxonomy"] == kraken_name]
    if not match.empty:
        summary["Detected_by_Krakenuniq"] = True
        summary["Krakenuniq_reads"] = match[reads_col].values[0]
        summary["Escore"] = match["Escore"].values[0]
        # Add Guellil et al Evalue if available
        if "Guellil_et_al_Evalue" in match.columns:
            summary["Guellil_et_al_Evalue"] = match["Guellil_et_al_Evalue"].values[0]

    # HOPS match
    if sample_col is not None and not hops_df.empty:
        match = hops_df[hops_df["Species"] == hops_name_safe]
        if not match.empty and sample_col in match.columns:
            hops_count = match[sample_col].values[0]
            summary["HOPS_score"] = hops_count
            if hops_count > 1:
                summary["Detected_by_HOPS"] = True

    # BAM path - try multiple locations (merged sample-level, per-PCR, legacy)
    bam_path = None
    possible_paths = [
        os.path.join(args.bam_dir, f"{args.sample}_{pathogen_safe}.dedup.merged.bam"),  # Sample-level merged
        os.path.join(args.bam_dir, f"{args.sample}_{pathogen_safe}.dedup.bam"),  # Per-PCR or sample-level
        f"results/pathogen/{args.sample}/pathogen_mapping/{args.sample}_{pathogen_safe}.dedup.bam",  # Direct path
        os.path.join(args.bam_dir, f"{args.sample}_{pathogen_safe}_F4_q30_sort.bam"),  # Legacy
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            bam_path = path
            break

    if bam_path and os.path.exists(bam_path):
        try:
            bam = pysam.AlignmentFile(bam_path, "rb")
            # Count mapped reads only
            mapped_count = sum(1 for read in bam.fetch(until_eof=True) if not read.is_unmapped)
            summary["Mapped_reads"] = mapped_count
            
            # Calculate read mapping ratio: Mapped_reads / Krakenuniq_reads
            # This metric indicates consistency between KrakenUniq detection and BWA/Bowtie2 mapping
            # > 1.0: More reads mapped than KrakenUniq detected (good - suggests KrakenUniq was conservative)
            # < 1.0: Fewer reads mapped than KrakenUniq detected (could indicate contamination or misclassification)
            # Close to 1.0: Perfect agreement between methods
            # Compute Read_mapping_ratio robustly (Krakenuniq_reads can be numpy types or strings)
            kraken_reads_raw = summary.get("Krakenuniq_reads", 0)
            try:
                kraken_reads_num = float(kraken_reads_raw)
            except Exception:
                kraken_reads_num = 0.0

            if kraken_reads_num > 0:
                mapping_ratio = mapped_count / kraken_reads_num
                summary["Read_mapping_ratio"] = round(mapping_ratio, 4)
            else:
                summary["Read_mapping_ratio"] = "NA"
            
            bam.close()

            # Get read lengths
            bam = pysam.AlignmentFile(bam_path, "rb")
            lengths = []
            for read in bam.fetch(until_eof=True):
                if not read.is_unmapped and read.query_length is not None and read.query_length > 0:
                    lengths.append(read.query_length)
            bam.close()

            if lengths:
                summary["Read_length_mean"] = round(np.mean(lengths), 2)
                summary["Read_length_median"] = round(np.median(lengths), 2)
            else:
                summary["Read_length_mean"] = "NA"
                summary["Read_length_median"] = "NA"
        except Exception as e:
            print(f"[WARNING] Error processing BAM {bam_path}: {e}")
            summary["Mapped_reads"] = "error"
            summary["Read_mapping_ratio"] = "error"
            summary["Read_length_mean"] = "error"
            summary["Read_length_median"] = "error"
    else:
        summary["Mapped_reads"] = "NA"
        summary["Read_mapping_ratio"] = "NA"
        summary["Read_length_mean"] = "NA"
        summary["Read_length_median"] = "NA"

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

    # Relative Entropy from files (pathopipe method: two window sizes)
    # Read both 100bp and 1000bp entropy values
    entropy_100_file = os.path.join(args.qualimap_dir, f"{args.sample}_{pathogen_safe}.entropy_100bp.txt")
    entropy_1000_file = os.path.join(args.qualimap_dir, f"{args.sample}_{pathogen_safe}.entropy_1000bp.txt")
    entropy_legacy_file = os.path.join(args.qualimap_dir, f"{args.sample}_{pathogen_safe}.mean_entropy.txt")
    
    entropy_100_value = None
    entropy_1000_value = None
    
    # Try to read 100bp entropy
    if os.path.exists(entropy_100_file):
        try:
            with open(entropy_100_file) as f:
                entropy_100_value = float(f.readline().strip())
                summary["Relative_entropy_100bp"] = round(entropy_100_value, 4)
        except Exception:
            summary["Relative_entropy_100bp"] = "error"
    else:
        summary["Relative_entropy_100bp"] = "NA"
    
    # Try to read 1000bp entropy
    if os.path.exists(entropy_1000_file):
        try:
            with open(entropy_1000_file) as f:
                entropy_1000_value = float(f.readline().strip())
                summary["Relative_entropy_1000bp"] = round(entropy_1000_value, 4)
        except Exception:
            summary["Relative_entropy_1000bp"] = "error"
    else:
        summary["Relative_entropy_1000bp"] = "NA"
    
    # For backward compatibility, also store as "Relative_entropy" (use 1000bp value)
    if entropy_1000_value is not None:
        summary["Relative_entropy"] = round(entropy_1000_value, 4)
    elif os.path.exists(entropy_legacy_file):
        # Fallback to legacy file if new files don't exist
        try:
            with open(entropy_legacy_file) as f:
                summary["Relative_entropy"] = float(f.readline().strip())
        except Exception:
            summary["Relative_entropy"] = "error"
    else:
        summary["Relative_entropy"] = "NA"
    
    # Determine if this is a virus for threshold selection
    # According to Nature paper: ALL viruses use 0.7 threshold (not just short ones)
    # Detection: check if pathogen name contains "virus" (case-insensitive)
    is_virus = False
    if "virus" in kraken_name.lower():
        is_virus = True
    # Also check BAM file to get reference length as additional check
    # (short genomes < 10kb are also treated as viruses)
    bam_file = os.path.join(args.bam_dir, f"{args.sample}_{pathogen_safe}.dedup.bam")
    if os.path.exists(bam_file):
        try:
            import pysam
            bam = pysam.AlignmentFile(bam_file, "rb")
            if len(bam.references) > 0:
                ref_len = bam.lengths[0]
                # Short genomes: < 10kb (also use virus threshold)
                if ref_len < 10000:
                    is_virus = True
            bam.close()
        except Exception:
            pass  # If we can't determine, rely on name-based detection

    # Breadth ratio (observed/expected) from MappingStats outputs (sample-level)
    # This is the ratio of observed to expected breadth of coverage (as in Nature paper)
    breadth_ratio_file = os.path.join(args.qualimap_dir, f"{args.sample}_{pathogen_safe}.breadth_ratio.txt")
    if os.path.exists(breadth_ratio_file):
        try:
            val = float(open(breadth_ratio_file).readline().strip())
            summary["Breadth_ratio"] = round(val, 4)
        except Exception:
            summary["Breadth_ratio"] = "error"

    # ANI from file (sample-level)
    ani_file = os.path.join(args.qualimap_dir, f"{args.sample}_{pathogen_safe}.ani.txt")
    ani_value = None
    if os.path.exists(ani_file):
        try:
            with open(ani_file) as f:
                for line in f:
                    if "ANI" in line:
                        ani_str = line.strip().split("≈")[-1].strip()
                        summary["ANI"] = ani_str
                        # Extract numeric value for scoring
                        try:
                            ani_value = float(ani_str.replace("%", "").strip())
                        except:
                            pass
                        break
        except Exception:
            summary["ANI"] = "error"

    # Edit distance metrics from file (sample-level)
    edit_r2_file = os.path.join(args.qualimap_dir, f"{args.sample}_{pathogen_safe}.edit_distance_logr2.txt")
    edit_r2_damage_file = os.path.join(
        args.qualimap_dir,
        f"{args.sample}_{pathogen_safe}.edit_distance_logr2_damaged.txt",
    )
    edit_r2_default_file = os.path.join(
        args.qualimap_dir,
        f"{args.sample}_{pathogen_safe}.edit_distance_logr2_default.txt",
    )
    # Legacy filenames (pre-renamed runs)
    edit_r2_default_legacy = os.path.join(
        args.qualimap_dir,
        f"{args.sample}_{pathogen_safe}.edit_distance_logr2_no_damage.txt",
    )

    decay_quality_score = None  # Decay Quality Score (damaged-read subset preferred)
    decay_quality_score_default = None
    edit_r2_value = None   # legacy 0–1 transformed metric (optional)

    def _first_line_decay_quality(path: str):
        with open(path) as f:
            for ln in f:
                ln = ln.strip()
                if ln == "":
                    continue
                if ln == "NA":
                    return None
                return float(ln)
        return None

    # Prefer damaged vs default (non-damaged) subset metrics when both files are present.
    default_metrics = edit_r2_default_file
    if not os.path.exists(default_metrics) and os.path.exists(edit_r2_default_legacy):
        default_metrics = edit_r2_default_legacy
    if os.path.exists(edit_r2_damage_file) and os.path.exists(default_metrics):
        try:
            decay_quality_score = _first_line_decay_quality(edit_r2_damage_file)
            decay_quality_score_default = _first_line_decay_quality(default_metrics)

            if decay_quality_score is not None:
                summary["Edit_distance_decay_quality"] = round(decay_quality_score, 4)
            else:
                summary["Edit_distance_decay_quality"] = "NA"

            if decay_quality_score_default is not None:
                summary["Edit_distance_decay_quality_default"] = round(
                    decay_quality_score_default, 4
                )
        except Exception:
            summary["Edit_distance_decay_quality"] = "error"
            summary["Edit_distance_decay_quality_default"] = "error"
    elif os.path.exists(edit_r2_file):
        try:
            with open(edit_r2_file) as f:
                lines = [ln.strip() for ln in f if ln.strip() != ""]
                # First line is now the Decay Quality Score (primary metric)
                if len(lines) >= 1 and lines[0] != "NA":
                    try:
                        decay_quality_score = float(lines[0])
                    except Exception:
                        decay_quality_score = None
                # Scan remaining lines for other metrics
                for ln in lines[1:]:
                    if ln.startswith("Legacy_metric="):
                        try:
                            edit_r2_value = float(ln.split("=", 1)[1])
                        except Exception:
                            edit_r2_value = None
                    elif ln.startswith("Decay_quality_score="):
                        # Also check explicit line (for robustness)
                        if decay_quality_score is None:
                            try:
                                decay_quality_score = float(ln.split("=", 1)[1])
                            except Exception:
                                pass

                # Store Decay Quality Score as the main summary metric
                if decay_quality_score is not None:
                    summary["Edit_distance_decay_quality"] = round(decay_quality_score, 4)
                else:
                    summary["Edit_distance_decay_quality"] = "NA"
        except Exception:
            summary["Edit_distance_decay_quality"] = "error"

    # Genus ranking from file (sample-level)
    genus_ranking_file = os.path.join(args.qualimap_dir, f"{args.sample}_{pathogen_safe}.genus_ranking.txt")
    genus_ranking_value = None
    if os.path.exists(genus_ranking_file):
        try:
            with open(genus_ranking_file) as f:
                val = f.readline().strip()
                if val != "NA":
                    genus_ranking_value = int(val)
                    summary["Genus_ranking"] = f"#{genus_ranking_value}"  # Format as "#5" etc.
                else:
                    summary["Genus_ranking"] = "NA"
        except Exception:
            summary["Genus_ranking"] = "error"

    # Calculate Detection Score
    # Criteria: E-Score threshold, Reads threshold, HOPS (3 criteria if enabled), ANI, Damage, Breadth ratio, Entropy
    detection_score = 0
    criteria_passed = {}
    
    # Get thresholds from spreadsheet (pathogen-specific)
    pathogen_row = spreadsheet[spreadsheet["Krakenuniq name"] == kraken_name]
    
    # Get E-score threshold from spreadsheet, fallback to config default if missing
    if not pathogen_row.empty and "min_escore" in pathogen_row.columns:
        min_escore_val = pathogen_row["min_escore"].values[0]
        # Handle NaN/NA values
        if pd.notna(min_escore_val):
            min_escore = float(min_escore_val)
        else:
            min_escore = criteria_config.get("escore_threshold", 5)  # Fallback to config default
    else:
        min_escore = criteria_config.get("escore_threshold", 5)  # Fallback to config default
    
    # Get Reads threshold from spreadsheet, fallback to config default if missing
    if not pathogen_row.empty and "min_reads" in pathogen_row.columns:
        min_reads_val = pathogen_row["min_reads"].values[0]
        # Handle NaN/NA values
        if pd.notna(min_reads_val):
            min_reads = int(min_reads_val)
        else:
            min_reads = criteria_config.get("reads_threshold", 50)  # Fallback to config default
    else:
        min_reads = criteria_config.get("reads_threshold", 50)  # Fallback to config default
    
    # 1. Primary KrakenUniq + detection criterion:
    #    Either E-score-based (default) or Guellil et al. E-value-based depending on config.
    use_evalue = bool(criteria_config.get("use_evalue_for_detection", False))
    
    # Get E-value threshold from spreadsheet (per-pathogen) or config default
    max_evalue = criteria_config.get("guellil_evalue_threshold", 0.001)
    if not pathogen_row.empty:
        # Allow several possible column names for flexibility
        for col in ["Guellil_et_al_Evalue_threshold", "max_evalue", "evalue_threshold"]:
            if col in pathogen_row.columns:
                val = pathogen_row[col].values[0]
                if pd.notna(val):
                    try:
                        max_evalue = float(val)
                    except Exception:
                        pass
                break
    
    if use_evalue:
        # E-value criterion (larger is better)
        evalue_val = summary.get("Guellil_et_al_Evalue", "NA")
        evalue_label = f"E-value > {max_evalue}"
        try:
            evalue_num = float(evalue_val) if evalue_val != "NA" and evalue_val != "error" else None
            if evalue_num is not None and evalue_num > max_evalue:
                detection_score += 1
                criteria_passed[evalue_label] = True
            else:
                criteria_passed[evalue_label] = False
        except Exception:
            criteria_passed[evalue_label] = False
    else:
        # E-score criterion (larger is better)
        escore_val = summary.get("Escore", "NA")
        escore_label = f"E-Score >= {min_escore}"
        try:
            escore_num = float(escore_val) if escore_val != "NA" and escore_val != "error" else None
            if escore_num is not None and escore_num >= min_escore:
                detection_score += 1
                criteria_passed[escore_label] = True
            else:
                criteria_passed[escore_label] = False
        except Exception:
            criteria_passed[escore_label] = False
    
    # 2. Reads threshold (1 point) - check against pathogen-specific threshold from spreadsheet
    reads_val = summary.get("Krakenuniq_reads", 0)
    reads_label = f"Reads >= {min_reads}"
    try:
        reads_num = int(reads_val) if reads_val != "NA" and reads_val != "error" else 0
        if reads_num > 0 and reads_num >= min_reads:
            detection_score += 1
            criteria_passed[reads_label] = True
        else:
            criteria_passed[reads_label] = False
    except:
        criteria_passed[reads_label] = False
    
    # 3-5. HOPS criteria (only if HOPS is enabled)
    hops_enabled = sample_col is not None and not hops_df.empty
    hops_detection_thresh = criteria_config.get("hops_detection_threshold", 2)
    hops_edit_thresh = criteria_config.get("hops_edit_distance_threshold", 3)
    hops_damage_thresh = criteria_config.get("hops_damage_threshold", 4)
    
    if hops_enabled:
        # HOPS detection (1 point)
        if summary["Detected_by_HOPS"] and summary["HOPS_score"] >= hops_detection_thresh:
            detection_score += 1
            criteria_passed["HOPS Detection"] = True
        else:
            criteria_passed["HOPS Detection"] = False
        
        # HOPS edit distance (1 point)
        if summary["HOPS_score"] >= hops_edit_thresh:
            detection_score += 1
            criteria_passed["HOPS Edit Distance"] = True
        else:
            criteria_passed["HOPS Edit Distance"] = False
        
        # HOPS damage (1 point)
        if summary["HOPS_score"] >= hops_damage_thresh:
            detection_score += 1
            criteria_passed["HOPS Damage"] = True
        else:
            criteria_passed["HOPS Damage"] = False
    else:
        criteria_passed["HOPS Detection"] = None  # Not applicable
        criteria_passed["HOPS Edit Distance"] = None
        criteria_passed["HOPS Damage"] = None
    
    # 6. ANI threshold (1 point) - from config
    ani_threshold = criteria_config.get("ani_threshold", 96.5)
    ani_label = f"ANI > {ani_threshold}%"
    if ani_value is not None and ani_value > ani_threshold:
        detection_score += 1
        criteria_passed[ani_label] = True
    else:
        criteria_passed[ani_label] = False
    
    # 7. 5' C>T damage threshold (1 point) - from config
    damage_5p_thresh = criteria_config.get("damage_5p_ct_threshold", 0.01)
    damage_label = f"5' C>T Damage >= {damage_5p_thresh}"
    try:
        damage_5p = float(summary["Damage_5p_CtoT"]) if summary["Damage_5p_CtoT"] != "NA" and summary["Damage_5p_CtoT"] != "error" else None
        if damage_5p is not None and damage_5p >= damage_5p_thresh:
            detection_score += 1
            criteria_passed[damage_label] = True
        else:
            criteria_passed[damage_label] = False
    except:
        criteria_passed[damage_label] = False
    
    # 8. Breadth ratio threshold (1 point) - from config
    breadth_thresh = criteria_config.get("breadth_ratio_threshold", 0.8)
    breadth_label = f"Breadth Ratio >= {breadth_thresh}"
    try:
        breadth_ratio = float(summary["Breadth_ratio"]) if summary["Breadth_ratio"] != "NA" and summary["Breadth_ratio"] != "error" else None
        if breadth_ratio is not None and breadth_ratio >= breadth_thresh:
            detection_score += 1
            criteria_passed[breadth_label] = True
        else:
            criteria_passed[breadth_label] = False
    except:
        criteria_passed[breadth_label] = False
    
    # 9. Relative entropy threshold (1 point) - pathopipe method
    # Use 1000bp window entropy for scoring (as in Nature paper)
    # Threshold: 0.9 for bacteria, 0.7 for viruses (all viruses, not just short ones)
    if is_virus:
        entropy_thresh = criteria_config.get("entropy_threshold_virus", 0.7)
        entropy_label = f"Entropy (1000bp) >= {entropy_thresh} (virus)"
    else:
        entropy_thresh = criteria_config.get("entropy_threshold", 0.9)
        entropy_label = f"Entropy (1000bp) >= {entropy_thresh} (bacteria)"
    
    try:
        # Use 1000bp entropy value for scoring
        entropy = entropy_1000_value
        if entropy is None:
            # Fallback to legacy Relative_entropy column
            entropy = float(summary["Relative_entropy"]) if summary["Relative_entropy"] != "NA" and summary["Relative_entropy"] != "error" else None
        
        if entropy is not None and entropy >= entropy_thresh:
            detection_score += 1
            criteria_passed[entropy_label] = True
        else:
            criteria_passed[entropy_label] = False
    except:
        criteria_passed[entropy_label] = False
    
    # 10. Edit distance criteria (damage vs no-damage when available)
    # Legacy behaviour: one metric "Edit Distance Decay Quality >= threshold" (+1 point)
    # Damage-split behaviour: two metrics (+0..2 points)
    #   - damage subset:     "Edit Distance Damage Decay Quality >= threshold" (+1)
    #   - no-damage subset:  "Edit Distance No-Damage Decay Quality >= min" (+1) [canonical]
    #                        or legacy: "<= edit_distance_no_damage_decay_quality_max" (+1)
    edit_decay_thresh = criteria_config.get("edit_distance_decay_quality_threshold", 0.65)
    no_damage_min = criteria_config.get("edit_distance_no_damage_decay_quality_min")
    no_damage_max = criteria_config.get("edit_distance_no_damage_decay_quality_max")
    edit_points_possible = 1  # legacy default

    try:
        default_val = decay_quality_score_default
        if decay_quality_score is not None and default_val is not None:
            edit_points_possible = 2
            damage_label = f"Edit distance damaged >= {edit_decay_thresh}"

            damage_pass = decay_quality_score >= edit_decay_thresh
            if no_damage_min is not None:
                ndm = float(no_damage_min)
                default_label = f"Edit distance default (non-damaged) >= {ndm}"
                default_pass = default_val >= ndm
            elif no_damage_max is not None:
                ndx = float(no_damage_max)
                default_label = f"Edit distance default (non-damaged) <= {ndx}"
                default_pass = default_val <= ndx
            else:
                ndm = 0.55
                default_label = f"Edit distance default (non-damaged) >= {ndm}"
                default_pass = default_val >= ndm

            if damage_pass:
                detection_score += 1
                criteria_passed[damage_label] = True
            else:
                criteria_passed[damage_label] = False

            if default_pass:
                detection_score += 1
                criteria_passed[default_label] = True
            else:
                criteria_passed[default_label] = False
        else:
            passed = (
                decay_quality_score is not None
                and decay_quality_score >= edit_decay_thresh
            )

            edit_r2_label = f"Edit Distance Decay Quality >= {edit_decay_thresh}"
            if passed:
                detection_score += 1
                criteria_passed[edit_r2_label] = True
            else:
                criteria_passed[edit_r2_label] = False
    except Exception:
        # If anything goes wrong, keep legacy label unset (or set to fail)
        edit_r2_label = f"Edit Distance Decay Quality >= {edit_decay_thresh}"
        criteria_passed[edit_r2_label] = False
    
    # 11. Genus ranking (1 point) - rank 1 is best (top position in genus)
    genus_ranking_label = "Genus Ranking = 1"
    try:
        if genus_ranking_value is not None and genus_ranking_value == 1:
            detection_score += 1
            criteria_passed[genus_ranking_label] = True
        else:
            criteria_passed[genus_ranking_label] = False
    except:
        criteria_passed[genus_ranking_label] = False
    
    # 12. Read mapping ratio threshold (1 point) - consistency between KrakenUniq and mapping
    read_ratio_thresh = criteria_config.get("read_mapping_ratio_threshold", 0.5)
    read_ratio_label = f"Read Mapping Ratio >= {read_ratio_thresh}"
    try:
        # Get read mapping ratio from summary (may be computed or from file)
        read_ratio_val = summary.get("Read_mapping_ratio", "NA")
        if read_ratio_val == "NA" or read_ratio_val == "error":
            # Try to compute from Mapped_reads and Krakenuniq_reads if ratio not available
            mapped_reads = summary.get("Mapped_reads", "NA")
            kraken_reads = summary.get("Krakenuniq_reads", 0)
            if mapped_reads != "NA" and mapped_reads != "error" and kraken_reads > 0:
                try:
                    read_ratio_val = float(mapped_reads) / float(kraken_reads)
                except:
                    read_ratio_val = None
            else:
                read_ratio_val = None
        else:
            read_ratio_val = float(read_ratio_val) if read_ratio_val != "NA" and read_ratio_val != "error" else None
        
        if read_ratio_val is not None and read_ratio_val >= read_ratio_thresh:
            detection_score += 1
            criteria_passed[read_ratio_label] = True
        else:
            criteria_passed[read_ratio_label] = False
    except:
        criteria_passed[read_ratio_label] = False
    
    # Calculate max possible score:
    # - base score assumes edit-distance contributes +1
    # - if damage/no-damage metrics are available, edit-distance contributes +2
    base_max = 9 if not hops_enabled else 12
    max_score = base_max + (edit_points_possible - 1)
    
    # Format score as "X/Y" (e.g., "3/6" or "5/9")
    summary["Score"] = f"{detection_score}/{max_score}"
    
    # Store only passed criteria as comma-separated list (cleaner than JSON)
    passed_criteria = [criterion for criterion, passed in criteria_passed.items() if passed is True]
    summary["Criteria_passed"] = ", ".join(passed_criteria) if passed_criteria else "None"

    rows.append(summary)

# -------------------- Save summary --------------------
os.makedirs(os.path.dirname(args.output), exist_ok=True)
if rows:
    df = pd.DataFrame(rows)
else:
    # Never write a 0-byte CSV: create a header-only file with the expected columns.
    df = pd.DataFrame(
        columns=[
            "Sample",
            "Pathogen",
            "Host_best_species",
            "Host_aligner",
            "Pathogen_screening_only",
            "Krakenuniq_reads",
            "Escore",
            "Guellil_et_al_Evalue",
            "Detected_by_Krakenuniq",
            "HOPS_score",
            "Detected_by_HOPS",
            "Mapped_reads",
            "Read_mapping_ratio",
            "Coverage",
            "Evenness",
            "Read_length_mean",
            "Read_length_median",
            "Damage_5p_CtoT",
            "Relative_entropy",
            "Relative_entropy_100bp",
            "Relative_entropy_1000bp",
            "Breadth_ratio",
            "ANI",
            "Edit_distance_decay_quality",
            "Edit_distance_decay_quality_default",
            "Genus_ranking",
            "Score",
            "Criteria_passed",
        ]
    )
df.to_csv(args.output, index=False)

# -------------------- Cleanup intermediate files if enabled --------------------
if args.cleanup:
    print(f"[CLEANUP] Removing intermediate pathogen files for sample {args.sample}...")
    # Only clean up files for pathogens that were actually detected (in observed set)
    for kraken_name in observed:
        pathogen_safe = safe_pathogen_name(kraken_name)
        
        # Remove intermediate pathogen BAMs (before dedup) - check both sample_pathogen_mapping and per-PCR locations
        intermediate_files = [
            # Per-sample merged intermediate files (if they exist)
            os.path.join(args.bam_dir, f"{args.sample}_{pathogen_safe}.sai"),
            # Per-PCR intermediate files (need to check all PCRs for this sample)
            # We'll use glob to find them
        ]
        
        # Use glob to find per-PCR intermediate files
        per_pcr_patterns = [
            f"results/pathogen/*/pathogen_mapping/*_{pathogen_safe}.sai",
            f"results/pathogen/*/pathogen_mapping/*_{pathogen_safe}_F4.bam",
            f"results/pathogen/*/pathogen_mapping/*_{pathogen_safe}_F4_q30.bam",
            f"results/pathogen/*/pathogen_mapping/*_{pathogen_safe}_F4_q30.sorted.bam",
            f"results/libraries/*/pathogen_mapping/*_{pathogen_safe}.sai",
            f"results/libraries/*/pathogen_mapping/*_{pathogen_safe}_F4.bam",
            f"results/libraries/*/pathogen_mapping/*_{pathogen_safe}_F4_q30.bam",
            f"results/libraries/*/pathogen_mapping/*_{pathogen_safe}_F4_q30.sorted.bam",
        ]
        
        for pattern in per_pcr_patterns:
            intermediate_files.extend(glob.glob(pattern))
        
        for filepath in intermediate_files:
            if os.path.exists(filepath):
                try:
                    os.remove(filepath)
                    print(f"  Removed: {filepath}")
                except Exception as e:
                    print(f"  Warning: Could not remove {filepath}: {e}")
    
    print(f"[CLEANUP] Intermediate pathogen file cleanup complete for {args.sample}.")
