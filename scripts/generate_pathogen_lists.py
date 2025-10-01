# scripts/generate_pathogen_lists.py
import pandas as pd
import sys
from pathlib import Path

spreadsheet = sys.argv[1]
df = pd.read_csv(spreadsheet)

# Create the output directory if it doesn't exist
Path("lists").mkdir(exist_ok=True)

# Map lowercase column names to real names
col_map = {col.lower().strip(): col for col in df.columns}

def match_column(possible_keys, available_keys):
    for key in available_keys:
        for target in possible_keys:
            if target in key:
                return available_keys[key]
    return None

kraken_col = match_column(["krakenuniq"], col_map)
hops_col = match_column(["hops"], col_map)

if kraken_col is None or hops_col is None:
    raise KeyError(
        f"❌ Could not find required columns in spreadsheet.\n"
        f"  ➤ Found columns: {list(df.columns)}\n"
        f"  ➤ Needed something like 'Krakenuniq name' and 'Hops name'"
    )

# Write the actual data to the correct paths
df[kraken_col].dropna().drop_duplicates().str.strip().to_csv(
    "lists/krakenuniq_pathogen_list.txt", index=False, header=False
)

df[hops_col].dropna().drop_duplicates().str.strip().to_csv(
    "lists/hops_pathogen_list.txt", index=False, header=False
)
