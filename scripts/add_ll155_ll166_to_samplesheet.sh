#!/bin/bash
# Script to add LL155-LL166 entries to samples_LL110-LL166.tsv
# Usage: ./add_ll155_ll166_to_samplesheet.sh

SAMPLE_SHEET="/raid_md0/louis/screening/samples_LL110-LL166.tsv"
BASE_DIR="/raid_md0/louis/screening/D-MCRG-Mar25-X10B-2/raw_data"
SEQ_RUN="D-MCRG-Mar25-X10B-2"

# Check if sample sheet exists
if [[ ! -f "$SAMPLE_SHEET" ]]; then
    echo "ERROR: Sample sheet not found: $SAMPLE_SHEET"
    exit 1
fi

# LL155-LL166 samples
samples=(
    "LL155-6-30"
    "LL156-6-31"
    "LL157-6-32"
    "LL158-6-33"
    "LL159-6-34"
    "LL160-6-35"
    "LL161-6-36"
    "LL162-6-37"
    "LL163-7-38"
    "LL164-7-39"
    "LL165-7-40"
    "LL166-7-41"
)

echo "Adding entries to $SAMPLE_SHEET..."

# Extract sample IDs (first part before dash) and add entries
for entry in "${samples[@]}"; do
    sample_id=$(echo "$entry" | cut -d'-' -f1)
    pcr_id="$entry"
    r1="${BASE_DIR}/${entry}_R1.fastq.gz"
    r2="${BASE_DIR}/${entry}_R2.fastq.gz"
    rglb="$pcr_id"
    
    # Verify files exist
    if [[ ! -f "$r1" ]]; then
        echo "WARNING: R1 file not found: $r1"
        continue
    fi
    if [[ ! -f "$r2" ]]; then
        echo "WARNING: R2 file not found: $r2"
        continue
    fi
    
    # Append to sample sheet (tab-separated, matching existing format)
    printf "%s\t%s\t%s\t%s\t%s\t%s\tLOCAL\n" \
        "$sample_id" "$pcr_id" "$r1" "$r2" "$rglb" "$SEQ_RUN" >> "$SAMPLE_SHEET"
    
    echo "  Added: $sample_id -> $pcr_id"
done

echo ""
echo "Done! Added ${#samples[@]} entries to $SAMPLE_SHEET"
echo "You can verify with: tail -12 $SAMPLE_SHEET"

