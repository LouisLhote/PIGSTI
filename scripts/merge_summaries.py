import pandas as pd
import sys
from pathlib import Path

summaries = snakemake.input
output = snakemake.output[0]

dfs = []

for f in summaries:
    try:
        # Skip if file is completely empty
        if Path(f).stat().st_size == 0:
            print(f"⚠️ Skipping empty file: {f}", file=sys.stderr)
            continue
        
        df = pd.read_csv(f)
        if df.empty:
            print(f"⚠️ Skipping header-only file: {f}", file=sys.stderr)
            continue

        dfs.append(df)

    except pd.errors.EmptyDataError:
        print(f"❌ Skipping unreadable (empty) file: {f}", file=sys.stderr)
        continue
    except Exception as e:
        print(f"❌ Unexpected error with {f}: {e}", file=sys.stderr)
        raise

if dfs:
    merged = pd.concat(dfs, ignore_index=True)
    merged.to_excel(output, index=False)
    print(f"✅ Merged {len(dfs)} files to {output}")
else:
    print("❌ No valid input files to merge. Aborting.", file=sys.stderr)
    sys.exit(1)
