import pandas as pd
import sys
from pathlib import Path

summaries = snakemake.input
output = snakemake.output[0]

dfs = []
header_only_cols = None

for f in summaries:
    try:
        # Skip if file is completely empty
        if Path(f).stat().st_size == 0:
            print(f"Warning: Skipping empty file: {f}", file=sys.stderr)
            continue

        df = pd.read_csv(f)
        if df.empty:
            print(f"Warning: Skipping header-only file: {f}", file=sys.stderr)
            if header_only_cols is None and len(df.columns) > 0:
                header_only_cols = list(df.columns)
            continue

        dfs.append(df)

    except pd.errors.EmptyDataError:
        print(f"Warning: Skipping unreadable (empty) file: {f}", file=sys.stderr)
        continue
    except Exception as e:
        print(f"Error with {f}: {e}", file=sys.stderr)
        raise

sheet = "pathogen_summary"

if dfs:
    merged = pd.concat(dfs, ignore_index=True)
else:
    # No pathogen hits in any sample — still write a valid empty workbook so
    # downstream rules (heatmap, comprehensive summary) can continue.
    cols = header_only_cols or ["Sample", "Pathogen", "Score"]
    merged = pd.DataFrame(columns=cols)
    print(
        f"Warning: No pathogen rows to merge; writing empty sheet with columns: {cols}",
        file=sys.stderr,
    )

if output.endswith(".xlsx") or output.endswith(".xls"):
    merged.to_excel(output, index=False, sheet_name=sheet)
else:
    merged.to_csv(output, index=False)

print(f"Merged {len(dfs)} file(s) ({len(merged)} rows) to {output} (sheet={sheet})")
