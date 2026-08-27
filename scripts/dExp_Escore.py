import sys
import yaml
import pandas as pd
import numpy as np
from io import StringIO

def double_exp(x, a=1.3, b=18):
    try:
        return a ** (b * x)
    except OverflowError:
        return float("inf")

def e_score_dexp(nb_kmer, nb_read, cov):
    return (nb_kmer / nb_read) * double_exp(cov) if nb_read > 0 else 0

e_score_dexp_vec = np.vectorize(e_score_dexp)

# -----------------------------
# Input arguments
# -----------------------------
kraken_report = sys.argv[1]
output_genus = sys.argv[2]
output_species = sys.argv[3]
output_pathogen = sys.argv[4]
config_file = sys.argv[5]

# -----------------------------
# Load config and spreadsheet
# -----------------------------
with open(config_file) as f:
    config = yaml.safe_load(f)

default_min_reads = int(config["defaults"]["min_reads"])
default_min_escore = float(config["defaults"]["min_escore"])
criteria_cfg = config.get("pathogen_detection_criteria", {}) if isinstance(config, dict) else {}
use_evalue_for_detection = bool(criteria_cfg.get("use_evalue_for_detection", True))
default_evalue_threshold = float(criteria_cfg.get("guellil_evalue_threshold", 0.001))

spreadsheet = pd.read_csv("config/Pathogen_spreadsheet.csv")
spreadsheet.columns = [col.strip().lower() for col in spreadsheet.columns]
spreadsheet = spreadsheet.rename(columns={"tax id": "taxid"})
spreadsheet["taxid"] = spreadsheet["taxid"].astype(int)

# -----------------------------
# Load Kraken report (robust to empty / malformed files)
# -----------------------------
with open(kraken_report) as f:
    lines = f.readlines()

# Find the header line that starts with '%'. If none exists (e.g. empty or
# malformed Kraken report), short-circuit by writing empty outputs.
header_idx = None
for i, line in enumerate(lines):
    if line.startswith("%"):
        header_idx = i
        break

if header_idx is None:
    # No usable table in the Kraken report: create empty CSVs so downstream
    # rules see "no candidates" instead of crashing.
    empty_cols = [
        "reads",
        "uniq_kmers",
        "taxID",
        "rank",
        "taxonomy",
        "cov",
        "dup",
        "taxReads",
        "Escore",
        "Guellil_et_al_Evalue",
        "taxid",
        "min_escore",
        "min_reads",
    ]
    empty_df = pd.DataFrame(columns=empty_cols)
    # Genus and species outputs: simply empty tables.
    empty_df.to_csv(output_genus, index=False)
    empty_df.to_csv(output_species, index=False)
    # Pathogen output: also empty (no taxa pass thresholds because none exist).
    empty_df.to_csv(output_pathogen, index=False)
    sys.exit(0)

header = lines[header_idx].strip().lstrip("%").strip().split("\t")
data_lines = lines[header_idx + 1:]

df = pd.read_csv(StringIO("".join(data_lines)), sep="\t", names=header, engine="python")
df = df.rename(columns={"%": "pct", "kmers": "uniq_kmers", "taxName": "taxonomy"})

required_columns = [
    "reads", "uniq_kmers", "taxID", "rank", "taxonomy",
    "cov", "dup", "taxReads"
]
if "pct" in df.columns:
    required_columns.append("pct")

missing = set(required_columns) - set(df.columns)
if missing:
    raise KeyError(f"Missing expected column(s) in Kraken report: {missing}")

df = df[required_columns]
df = df.astype({
    "reads": int,
    "uniq_kmers": int,
    "taxID": int,
    "rank": str,
    "taxonomy": str,
    "cov": float,
    "dup": float,
    "taxReads": int,
    **({"pct": float} if "pct" in df.columns else {})
})

# -----------------------------
# Compute Escore
# -----------------------------
df["Escore"] = e_score_dexp_vec(df["uniq_kmers"], df["reads"], df["cov"])

# -----------------------------
# Compute Guellil et al Evalue: (Kmers/reads)*coverage
# -----------------------------
df["Guellil_et_al_Evalue"] = (df["uniq_kmers"] / df["reads"]) * df["cov"]
df["Guellil_et_al_Evalue"] = df["Guellil_et_al_Evalue"].replace([np.inf, -np.inf], np.nan)
df["Guellil_et_al_Evalue"] = df["Guellil_et_al_Evalue"].fillna(0)

# -----------------------------
# Write outputs (genus/species)
# -----------------------------
df.to_csv(output_genus, index=False)
df[df["rank"] == "species"].to_csv(output_species, index=False)

# -----------------------------
# Pathogen filtering (use only spreadsheet taxid)
# -----------------------------
known_taxids = spreadsheet["taxid"].unique()
# Merge thresholds from spreadsheet to ensure *only* those taxon are selected and proper fields used
threshold_cols = ["taxid", "min_escore", "min_reads"]
candidate_evalue_cols = [
    "guellil_et_al_evalue_threshold",
    "max_evalue",
    "evalue_threshold",
]
present_evalue_cols = [c for c in candidate_evalue_cols if c in spreadsheet.columns]
threshold_cols.extend(present_evalue_cols)

thresholds = spreadsheet[threshold_cols].drop_duplicates().copy()
thresholds = thresholds.fillna({
    "min_escore": default_min_escore,
    "min_reads": default_min_reads,
})

# Compute a single effective evalue threshold from allowed spreadsheet aliases.
if present_evalue_cols:
    thresholds["guellil_et_al_evalue_threshold"] = (
        thresholds[present_evalue_cols]
        .apply(pd.to_numeric, errors="coerce")
        .bfill(axis=1)
        .iloc[:, 0]
    )
else:
    thresholds["guellil_et_al_evalue_threshold"] = np.nan
thresholds["guellil_et_al_evalue_threshold"] = thresholds["guellil_et_al_evalue_threshold"].fillna(
    default_evalue_threshold
)

# Merge df (Kraken report) species with spreadsheet, keep only pathogens from spreadsheet
pathogen_df = df[df["taxID"].isin(known_taxids)].copy()
pathogen_df = pathogen_df.merge(thresholds, left_on="taxID", right_on="taxid", how="inner")

# The only "taxid" we use is from the spreadsheet (merge ensures this)
# Apply filters according to spreadsheet thresholds
pathogen_df["min_escore"] = pathogen_df["min_escore"].fillna(default_min_escore)
pathogen_df["min_reads"] = pathogen_df["min_reads"].fillna(default_min_reads)
pathogen_df["guellil_et_al_evalue_threshold"] = pathogen_df["guellil_et_al_evalue_threshold"].fillna(
    default_evalue_threshold
)

if use_evalue_for_detection:
    # Guellil metric mode (project-specific): larger values are better.
    filtered_pathogens = pathogen_df[
        (pathogen_df["reads"] >= pathogen_df["min_reads"]) &
        (pathogen_df["Guellil_et_al_Evalue"] > pathogen_df["guellil_et_al_evalue_threshold"])
    ].copy()
else:
    filtered_pathogens = pathogen_df[
        (pathogen_df["reads"] >= pathogen_df["min_reads"]) &
        (pathogen_df["Escore"] >= pathogen_df["min_escore"])
    ].copy()

# Choose and reorder final output columns—all downstream code expects lowercase 'taxid'
out_cols = [
    'reads', 'uniq_kmers', 'taxid', 'rank', 'taxonomy', 'cov', 'dup',
    'taxReads', 'Escore', 'Guellil_et_al_Evalue', 'min_escore', 'min_reads',
    'guellil_et_al_evalue_threshold'
]
out_cols = [c for c in out_cols if c in filtered_pathogens.columns]
filtered_pathogens[out_cols].to_csv(output_pathogen, index=False)
