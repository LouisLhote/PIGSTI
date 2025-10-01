# scripts/prepare_bwa_targets.py

import pandas as pd

# Get inputs/outputs directly from Snakemake API
escore_file = snakemake.input.escore
spreadsheet_file = snakemake.input.spreadsheet
output_file = snakemake.output.targets
sample = snakemake.wildcards.sample

# Read input CSVs
escore = pd.read_csv(escore_file)
spreadsheet = pd.read_csv(spreadsheet_file)

# Normalize column names
escore.columns = escore.columns.str.strip().str.lower()
spreadsheet.columns = spreadsheet.columns.str.strip().str.lower()

# Convert taxid to str (required for merging)
escore['taxid'] = escore['taxid'].astype(str)
spreadsheet['taxid'] = spreadsheet['taxid'].astype(str)

# Merge escore results with pathogen spreadsheet based on taxid
merged = pd.merge(escore, spreadsheet, on='taxid', how='inner')

# If no match found after merge, warn and write empty file
if merged.empty:
    print(f"[WARNING] No pathogens matched for sample {sample}. Writing empty output.")
    pd.DataFrame(columns=['taxid', 'Krakenuniq name', 'bwa index', 'ref_name']).to_csv(output_file, sep='\t', index=False)
else:
    # Sanitize ref_name for file-safe output paths
    merged['ref_name'] = merged['krakenuniq name'].str.replace(r'[^\w]+', '_', regex=True)

    # Output required columns
    merged[['taxid', 'krakenuniq name', 'bwa index', 'ref_name']].to_csv(output_file, sep='\t', index=False)
