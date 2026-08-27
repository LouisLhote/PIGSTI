#!/usr/bin/env bash
#
# Wrapper to:
#   1) Prepare and download ENA FASTQs using ena_to_pigsti_samplesheet.py
#   2) Merge ENA + local samples into config/samples.tsv
#   3) Run the full PIGSTI Snakemake pipeline
#
# Usage (default paths work with your current layout):
#   bash run_pigsti_with_ena.sh [snakemake-args...]
#
# You can customize:
#   - LOCAL_SAMPLES_TSV: path to your local-only samples (if any)
#   - ENA_FILEREPORT_GLOB: where ENA filereport TSVs live
#   - ENA_OUTPUT_DIR: where to download ENA FASTQs
#   - MERGED_SAMPLES_TSV: final samples.tsv PIGSTI will use

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"

# ---- User-configurable paths ----

# Local samples (your own FASTQs). If this file doesn't exist, only ENA samples will be used.
LOCAL_SAMPLES_TSV="${LOCAL_SAMPLES_TSV:-config/local_samples.tsv}"

# Where ENA filereport TSVs are stored.
# All matching files will be passed to ena_to_pigsti_samplesheet.py.
ENA_FILEREPORT_GLOB="${ENA_FILEREPORT_GLOB:-config/ENA/filereport_read_run_*.tsv}"

# Where to download ENA FASTQs
ENA_OUTPUT_DIR="${ENA_OUTPUT_DIR:-/raid_md0/louis/ENA_raw_pigsti}"

# Final samples.tsv used by PIGSTI
MERGED_SAMPLES_TSV="${MERGED_SAMPLES_TSV:-config/samples.tsv}"

# Temporary ENA-only samples file
ENA_SAMPLES_TSV="config/ena_samples.tsv"

# ---- Step 1: Find ENA filereport TSVs ----

shopt -s nullglob
ENA_FILES=( $ENA_FILEREPORT_GLOB )
shopt -u nullglob

if (( ${#ENA_FILES[@]} == 0 )); then
  echo "ERROR: No ENA filereport TSVs found matching: ${ENA_FILEREPORT_GLOB}" >&2
  echo "       Please place files like 'filereport_read_run_*.tsv' under config/ENA/ or adjust ENA_FILEREPORT_GLOB." >&2
  exit 1
fi

echo "Found ${#ENA_FILES[@]} ENA filereport(s):"
printf '  - %s\n' "${ENA_FILES[@]}"
echo

# ---- Step 2: Run ENA downloader + samples.tsv generator ----

mkdir -p "$(dirname "$ENA_SAMPLES_TSV")"
mkdir -p "$ENA_OUTPUT_DIR"

echo "=== [1/3] Generating ENA sample sheet and downloading FASTQs ==="
echo "Output ENA samples TSV: $ENA_SAMPLES_TSV"
echo "ENA FASTQ output dir:   $ENA_OUTPUT_DIR"
echo

python scripts/ena_to_pigsti_samplesheet.py \
  "${ENA_FILES[@]}" \
  --output-dir "$ENA_OUTPUT_DIR" \
  --samples-tsv "$ENA_SAMPLES_TSV"

echo
echo "ENA sample sheet created: $ENA_SAMPLES_TSV"
echo

if [[ ! -s "$ENA_SAMPLES_TSV" ]]; then
  echo "ERROR: ENA samples TSV is empty or missing: $ENA_SAMPLES_TSV" >&2
  exit 1
fi

# ---- Step 3: Merge local + ENA into config/samples.tsv with proper 'source' column ----

echo "=== [2/3] Merging local + ENA samples into $MERGED_SAMPLES_TSV ==="

mkdir -p "$(dirname "$MERGED_SAMPLES_TSV")"

ENA_HEADER="$(head -n1 "$ENA_SAMPLES_TSV")"

LOCAL_HEADER=""
if [[ -f "$LOCAL_SAMPLES_TSV" ]]; then
  LOCAL_HEADER="$(head -n1 "$LOCAL_SAMPLES_TSV")"
fi

if [[ -n "$LOCAL_HEADER" ]]; then
  # Determine final header: ensure it has a 'source' column
  if [[ "$LOCAL_HEADER" == *$'\tsource'* ]]; then
    HEADER="$LOCAL_HEADER"
    LOCAL_HAS_SOURCE=1
  else
    HEADER="${LOCAL_HEADER}\tsource"
    LOCAL_HAS_SOURCE=0
  fi
else
  # No local file; take header from ENA (which already includes 'source')
  HEADER="$ENA_HEADER"
  LOCAL_HAS_SOURCE=1
fi

printf '%s\n' "$HEADER" > "$MERGED_SAMPLES_TSV"

# Append local rows (if any)
if [[ -f "$LOCAL_SAMPLES_TSV" ]]; then
  echo "  - Adding local samples from $LOCAL_SAMPLES_TSV"
  if (( LOCAL_HAS_SOURCE == 1 )); then
    # Local already has a source column; copy as-is (without header)
    tail -n +2 "$LOCAL_SAMPLES_TSV" >> "$MERGED_SAMPLES_TSV"
  else
    # Add 'LOCAL' as source for each data row
    tail -n +2 "$LOCAL_SAMPLES_TSV" | awk -v OFS='\t' '{print $0, "LOCAL"}' >> "$MERGED_SAMPLES_TSV"
  fi
else
  echo "  - No local samples file found at $LOCAL_SAMPLES_TSV (only ENA samples will be used)."
fi

# Append ENA rows (skip header, ena_to_pigsti_samplesheet.py already sets source=ENA)
echo "  - Adding ENA samples from $ENA_SAMPLES_TSV"
tail -n +2 "$ENA_SAMPLES_TSV" >> "$MERGED_SAMPLES_TSV"

echo "Merged samples.tsv written to: $MERGED_SAMPLES_TSV"
echo

# ---- Step 4: Run the full PIGSTI pipeline ----

echo "=== [3/3] Running PIGSTI Snakemake pipeline ==="
echo "Using samples file: $MERGED_SAMPLES_TSV"
echo

snakemake --use-conda --cores 42 "$@"

echo
echo "PIGSTI run complete."


