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

spreadsheet = pd.read_csv("config/Pathogen_spreadsheet.csv")
spreadsheet.columns = [col.strip().lower() for col in spreadsheet.columns]
spreadsheet = spreadsheet.rename(columns={"tax id": "taxid"})
spreadsheet["taxid"] = spreadsheet["taxid"].astype(int)

# -----------------------------
# Load Kraken report
# -----------------------------
with open(kraken_report) as f:
    lines = f.readlines()

header_idx = next(i for i, line in enumerate(lines) if line.startswith("%"))
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
# Write outputs (genus/species)
# -----------------------------
df.to_csv(output_genus, index=False)
df[df["rank"] == "species"].to_csv(output_species, index=False)

# -----------------------------
# Pathogen filtering (use only spreadsheet taxid)
# -----------------------------
known_taxids = spreadsheet["taxid"].unique()
# Merge thresholds from spreadsheet to ensure *only* those taxon are selected and proper fields used
thresholds = spreadsheet[["taxid", "min_escore", "min_reads"]].drop_duplicates()
thresholds = thresholds.fillna({
    "min_escore": default_min_escore,
    "min_reads": default_min_reads
})

# Merge df (Kraken report) species with spreadsheet, keep only pathogens from spreadsheet
pathogen_df = df[df["taxID"].isin(known_taxids)].copy()
pathogen_df = pathogen_df.merge(thresholds, left_on="taxID", right_on="taxid", how="inner")

# The only "taxid" we use is from the spreadsheet (merge ensures this)
# Apply filters according to spreadsheet thresholds
pathogen_df["min_escore"] = pathogen_df["min_escore"].fillna(default_min_escore)
pathogen_df["min_reads"] = pathogen_df["min_reads"].fillna(default_min_reads)

filtered_pathogens = pathogen_df[
    (pathogen_df["reads"] >= pathogen_df["min_reads"]) &
    (pathogen_df["Escore"] >= pathogen_df["min_escore"])
].copy()

# Choose and reorder final output columns—all downstream code expects lowercase 'taxid'
out_cols = [
    'reads', 'uniq_kmers', 'taxid', 'rank', 'taxonomy', 'cov', 'dup',
    'taxReads', 'Escore', 'min_escore', 'min_reads'
]
out_cols = [c for c in out_cols if c in filtered_pathogens.columns]
filtered_pathogens[out_cols].to_csv(output_pathogen, index=False)
