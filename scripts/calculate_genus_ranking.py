#!/usr/bin/env python
"""
Parse KrakenUniq report to find the ranking position of a pathogen species
within its own genus. The ranking is based on read counts (descending order).
"""

import pandas as pd
import sys
import os
from io import StringIO

# Get inputs from Snakemake
try:
    kraken_report = snakemake.input.kraken_report
    spreadsheet = snakemake.input.spreadsheet
    pathogen_name = snakemake.params.pathogen_name
    output_ranking = snakemake.output.ranking
except AttributeError:
    # Fallback for command-line usage
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--kraken-report", required=True)
    parser.add_argument("--spreadsheet", required=True)
    parser.add_argument("--pathogen-name", required=True)
    parser.add_argument("--output-ranking", required=True)
    args = parser.parse_args()
    kraken_report = args.kraken_report
    spreadsheet = args.spreadsheet
    pathogen_name = args.pathogen_name
    output_ranking = args.output_ranking

# Load spreadsheet to get pathogen taxID
spreadsheet_df = pd.read_csv(spreadsheet)
spreadsheet_df["Krakenuniq name"] = spreadsheet_df["Krakenuniq name"].str.strip()

# Debug: print available pathogen names
sys.stderr.write(f"[DEBUG] Looking for pathogen: '{pathogen_name}'\n")
sys.stderr.write(f"[DEBUG] Available pathogens in spreadsheet: {list(spreadsheet_df['Krakenuniq name'].head(10))}\n")

# Find pathogen in spreadsheet
pathogen_row = spreadsheet_df[spreadsheet_df["Krakenuniq name"].str.lower() == pathogen_name.lower()]
if pathogen_row.empty:
    # Try exact match
    pathogen_row = spreadsheet_df[spreadsheet_df["Krakenuniq name"] == pathogen_name]

if pathogen_row.empty:
    sys.stderr.write(f"WARNING: Pathogen '{pathogen_name}' not found in spreadsheet. Available names: {list(spreadsheet_df['Krakenuniq name'].values[:5])}\n")
    with open(output_ranking, "w") as f:
        f.write("NA\n")
    sys.exit(0)

sys.stderr.write(f"[DEBUG] Found pathogen in spreadsheet: {pathogen_row['Krakenuniq name'].values[0]}\n")
sys.stderr.write(f"[DEBUG] Available columns in spreadsheet: {list(pathogen_row.columns)}\n")

# Get taxID - try different column name variations (case-insensitive)
taxid_col = None
possible_names = ["taxID", "Tax ID", "TaxID", "taxid", "tax_id", "Tax_id", "NCBI_taxid", "ncbi_taxid"]
# First try exact match
for col in possible_names:
    if col in pathogen_row.columns:
        taxid_col = col
        break
# If not found, try case-insensitive match
if taxid_col is None:
    col_lower_map = {col.lower(): col for col in pathogen_row.columns}
    for name in possible_names:
        if name.lower() in col_lower_map:
            taxid_col = col_lower_map[name.lower()]
            break

if taxid_col is None:
    sys.stderr.write(f"ERROR: Could not find Tax ID column in spreadsheet. Available columns: {list(pathogen_row.columns)}\n")
    with open(output_ranking, "w") as f:
        f.write("NA\n")
    sys.exit(1)

pathogen_taxid = int(pathogen_row[taxid_col].values[0])
pathogen_genus = pathogen_row["Genus"].values[0] if "Genus" in pathogen_row.columns else None

sys.stderr.write(f"[DEBUG] Pathogen taxID: {pathogen_taxid}, Genus: {pathogen_genus}\n")

# Load KrakenUniq report
if not os.path.exists(kraken_report):
    sys.stderr.write(f"ERROR: KrakenUniq report not found: {kraken_report}\n")
    with open(output_ranking, "w") as f:
        f.write("NA\n")
    sys.exit(1)

with open(kraken_report) as f:
    lines = f.readlines()

# Find header line
header_idx = None
for i, line in enumerate(lines):
    if line.startswith("%"):
        header_idx = i
        break

if header_idx is None:
    sys.stderr.write(f"ERROR: Could not find header in KrakenUniq report: {kraken_report}\n")
    with open(output_ranking, "w") as f:
        f.write("NA\n")
    sys.exit(1)

# Parse header
header = lines[header_idx].strip().lstrip("%").strip().split("\t")
data_lines = lines[header_idx + 1:]

# Read data
df = pd.read_csv(StringIO("".join(data_lines)), sep="\t", names=header, engine="python")
df = df.rename(columns={"taxName": "taxonomy"})

# Normalize rank values (handle single-letter codes like 'S', 'G', etc.)
rank_map = {
    "S": "species",
    "G": "genus",
    "F": "family",
    "O": "order",
    "C": "class",
    "P": "phylum",
    "K": "kingdom"
}
df["rank"] = df["rank"].astype(str).str.strip().map(lambda r: rank_map.get(r, r))

# Filter to species level, genus level, and family level (for virus fallback)
species_df = df[df["rank"] == "species"].copy()
genus_df = df[df["rank"] == "genus"].copy()
family_df = df[df["rank"] == "family"].copy()

# Find the genus taxID for our pathogen
# First, try to find the pathogen species in the report
pathogen_species_row = species_df[species_df["taxID"] == pathogen_taxid]

if pathogen_species_row.empty:
    # Some virus species may not be labeled strictly as 'species' in rank column
    # Fallback: look for any row with matching taxID regardless of rank
    pathogen_species_row = df[df["taxID"] == pathogen_taxid]

if pathogen_species_row.empty:
    # Pathogen not found in report
    sys.stderr.write(f"WARNING: Pathogen taxID {pathogen_taxid} not found in KrakenUniq report.\n")
    sys.stderr.write(f"[DEBUG] Available taxIDs in report (first 10): {list(df['taxID'].head(10).values)}\n")
    with open(output_ranking, "w") as f:
        f.write("NA\n")
    sys.exit(0)

sys.stderr.write(f"[DEBUG] Found pathogen in KrakenUniq report: {pathogen_species_row['taxonomy'].values[0]}\n")

# Get the genus taxID from the pathogen's taxonomy hierarchy
# In KrakenUniq reports, we need to find the parent genus of this species
pathogen_taxonomy = pathogen_species_row["taxonomy"].values[0].strip()
pathogen_taxid = int(pathogen_species_row["taxID"].values[0])

# Method 1: Look for the genus in the taxonomy hierarchy
# In KrakenUniq reports, species entries are indented under their genus
# We need to find the genus entry that appears before this species in the report
genus_taxid = None

# Find the integer position of the pathogen species in the full dataframe
# First try to find by rank "species", but also accept any row with matching taxID
# (viruses might have different rank labels)
pathogen_idx = None
for i, (idx, row) in enumerate(df.iterrows()):
    if int(row["taxID"]) == pathogen_taxid:
        # Found the pathogen - use this position for hierarchy lookup
        pathogen_idx = i  # Use enumerate position, not DataFrame index
        sys.stderr.write(f"[DEBUG] Found pathogen at position {i} with rank '{row['rank']}' and taxID {pathogen_taxid}\n")
        break

if pathogen_idx is not None:
    # Look backwards from the pathogen species to find the most recent genus entry
    # In KrakenUniq reports, species are listed under their parent genus
    # The hierarchy structure means the most recent genus before this species is its parent
    # For both bacteria and viruses, the hierarchy is: ... -> family -> genus -> species
    # So we should find the genus by scanning backwards and taking the first genus we encounter
    family_taxid = None  # Store for debugging, but viruses also have genus
    
    # Scan backwards from pathogen position to find parent genus
    # In KrakenUniq reports, the hierarchy is preserved in the order of entries
    # Children appear after their parents, so scanning backwards finds the parent
    for i in range(int(pathogen_idx) - 1, -1, -1):
        row = df.iloc[i]
        rank = str(row["rank"]).strip()
        
        # Normalize rank (handle single-letter codes like "G" for genus)
        rank_normalized = rank
        if rank in ["G", "g"]:
            rank_normalized = "genus"
        elif rank in ["F", "f"]:
            rank_normalized = "family"
        elif rank in ["S", "s"]:
            rank_normalized = "species"
        
        # If we find a genus, use it immediately (don't be too strict about name matching)
        # In KrakenUniq, the hierarchy structure ensures species are under their correct genus
        # This works for both bacteria and viruses
        if rank_normalized == "genus" or rank == "genus":
            genus_taxid = int(row["taxID"])
            genus_taxonomy = str(row["taxonomy"]).strip()
            sys.stderr.write(f"[DEBUG] Found genus taxID {genus_taxid} ({genus_taxonomy}) for pathogen '{pathogen_name}' (method: hierarchy lookup, rank='{rank}')\n")
            break
        elif (rank_normalized == "family" or rank == "family") and family_taxid is None:
            # Store family for debugging, but continue looking for genus
            # (Viruses also have genus, so we shouldn't stop here)
            family_taxid = int(row["taxID"])
            family_taxonomy = str(row["taxonomy"]).strip()
            sys.stderr.write(f"[DEBUG] Found family taxID {family_taxid} ({family_taxonomy}), continuing to look for genus...\n")
        elif rank in ["order", "class", "phylum", "kingdom", "superkingdom", "O", "C", "P", "K"]:
            # We've gone too far up the hierarchy, stop looking
            # But only stop if we haven't found a genus yet
            if genus_taxid is None:
                sys.stderr.write(f"[DEBUG] Reached {rank} level without finding genus, stopping search\n")
            break

# Method 2: If still not found, try to extract genus name from taxonomy string
# For viruses like "Cowpox virus", there might not be a traditional genus
# Try to find genus by matching the first word of the pathogen name
if genus_taxid is None:
    # Extract first word from pathogen name (e.g., "Cowpox" from "Cowpox virus")
    taxonomy_parts = pathogen_taxonomy.split()
    if len(taxonomy_parts) >= 1:
        potential_genus = taxonomy_parts[0].strip()
        # For viruses, the "genus" might be the virus name itself or a higher taxon
        # Try to find any genus entry that contains this word
        for idx, genus_row in genus_df.iterrows():
            genus_name = str(genus_row["taxonomy"]).strip()
            genus_first_word = genus_name.split()[0] if genus_name else ""
            # Match if first words match (case-insensitive)
            if potential_genus.lower() == genus_first_word.lower():
                genus_taxid = int(genus_row["taxID"])
                sys.stderr.write(f"[DEBUG] Found genus taxID {genus_taxid} by name matching: '{potential_genus}' -> '{genus_name}' (method: taxonomy string parsing)\n")
                break
        # If still not found, try substring match (for viruses that might be in a different format)
        if genus_taxid is None:
            for idx, genus_row in genus_df.iterrows():
                genus_name = str(genus_row["taxonomy"]).strip()
                # Check if potential_genus appears in genus name or vice versa
                if (potential_genus.lower() in genus_name.lower() or 
                    genus_name.lower() in potential_genus.lower()):
                    genus_taxid = int(genus_row["taxID"])
                    sys.stderr.write(f"[DEBUG] Found genus taxID {genus_taxid} by substring matching: '{potential_genus}' -> '{genus_name}' (method: taxonomy string parsing)\n")
                    break

# Method 3: If still not found and we have genus from spreadsheet, try that
if genus_taxid is None and pathogen_genus:
    genus_rows = genus_df[genus_df["taxonomy"].str.strip().str.lower() == pathogen_genus.lower()]
    if not genus_rows.empty:
        genus_taxid = int(genus_rows["taxID"].values[0])
        sys.stderr.write(f"[DEBUG] Found genus taxID {genus_taxid} from spreadsheet genus: '{pathogen_genus}'\n")

# Note: We no longer use family as fallback since viruses also have genus in the hierarchy
# If genus is still not found at this point, it means the hierarchy lookup failed
# and we should try the other methods (name matching, spreadsheet, etc.)
is_using_family_as_genus = False

if genus_taxid is None:
    sys.stderr.write(f"WARNING: Could not determine genus taxID for pathogen '{pathogen_name}'.\n")
    sys.stderr.write(f"[DEBUG] Pathogen taxonomy: '{pathogen_taxonomy}'\n")
    sys.stderr.write(f"[DEBUG] Available genera in report (first 10): {list(genus_df['taxonomy'].head(10).values)}\n")
    with open(output_ranking, "w") as f:
        f.write("NA\n")
    sys.exit(0)

# Find all species in the same genus
# In KrakenUniq reports, children follow their parent in the hierarchy
# For viruses, species taxonomy strings may NOT contain the genus name
# (e.g., "Sheeppox virus" doesn't contain "Capripoxvirus")
# So we collect all species that appear after the genus entry until the next genus
genus_species = []

# Find the position of the genus in the dataframe
genus_idx = None
for i, (idx, row) in enumerate(df.iterrows()):
    rank_check = str(row["rank"]).strip()
    is_genus = (rank_check == "genus" or rank_check in ["G", "g"])
    if is_genus and int(row["taxID"]) == genus_taxid:
        genus_idx = i
        break

if genus_idx is not None:
    # Collect all species that appear after this genus until we hit the next genus
    # This works for both bacteria and viruses
    for i in range(genus_idx + 1, len(df)):
        row = df.iloc[i]
        rank_check = str(row["rank"]).strip()
        is_genus = (rank_check == "genus" or rank_check in ["G", "g"])
        is_species = (rank_check == "species" or rank_check in ["S", "s"])
        
        if is_species:
            # Collect this species (it belongs to our genus by hierarchy position)
            taxonomy_str = row["taxonomy"].strip()
            genus_species.append({
                "taxID": int(row["taxID"]),
                "taxonomy": taxonomy_str,
                "reads": int(row["reads"])
            })
        elif is_genus:
            # We've moved to the next genus, stop collecting
            break
        elif rank_check in ["family", "F", "f", "order", "O", "o", "class", "C", "c"]:
            # We've moved up the hierarchy, stop collecting
            # (This shouldn't happen if hierarchy is correct, but safety check)
            break

# Fallback: If we didn't find species by hierarchy position, try name matching
# (This might work for bacteria where species names often contain genus)
if len(genus_species) == 0:
    genus_name = genus_df[genus_df["taxID"] == genus_taxid]["taxonomy"].values[0].strip()
    for idx, row in species_df.iterrows():
        taxonomy_str = row["taxonomy"].strip()
        # Check if taxonomy contains genus name (case-insensitive)
        # This is a fallback for cases where hierarchy position method didn't work
        if genus_name.lower() in taxonomy_str.lower():
            genus_species.append({
                "taxID": int(row["taxID"]),
                "taxonomy": taxonomy_str,
                "reads": int(row["reads"])
            })

# Sort by read count (descending)
genus_species.sort(key=lambda x: x["reads"], reverse=True)

# Find ranking position of our pathogen
ranking = None
for i, species in enumerate(genus_species, start=1):
    if species["taxID"] == pathogen_taxid:
        ranking = i
        break

if ranking is None:
    # Pathogen not found in genus species list
    sys.stderr.write(f"WARNING: Pathogen taxID {pathogen_taxid} not found in genus species list. Writing 'NA'.\n")
    with open(output_ranking, "w") as f:
        f.write("NA\n")
    sys.exit(0)

# Write ranking and all species in genus (for plotting)
with open(output_ranking, "w") as f:
    f.write(f"{ranking}\n")
    # Write all species in genus (one per line: taxID|taxonomy|reads)
    for species in genus_species:
        f.write(f"{species['taxID']}|{species['taxonomy']}|{species['reads']}\n")

print(f"✅ Found genus ranking: {ranking} for pathogen '{pathogen_name}' (taxID: {pathogen_taxid}) in genus (taxID: {genus_taxid})")
print(f"   Total species in genus: {len(genus_species)}")

