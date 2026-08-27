#!/bin/bash
# Script to remove duplicate entries from samples_LL110-LL166.tsv
# Keeps the first occurrence of each PCR ID

SAMPLE_SHEET="/raid_md0/louis/screening/samples_LL110-LL166.tsv"
BACKUP="${SAMPLE_SHEET}.backup_$(date +%Y%m%d_%H%M%S)"

# Create backup
echo "Creating backup: $BACKUP"
cp "$SAMPLE_SHEET" "$BACKUP"

# Remove duplicates using awk (keeps first occurrence)
echo "Removing duplicate entries..."
awk -F'\t' '
    BEGIN { OFS="\t" }
    NR == 1 { print; next }  # Print header
    {
        # Use PCR ID (column 2) as key, or sample (column 1) if PCR column is missing
        key = $2
        if (key == "" || key == "pcr") key = $1
        
        if (!seen[key]) {
            seen[key] = 1
            print
        } else {
            print "Skipping duplicate: " key > "/dev/stderr"
        }
    }
' "$SAMPLE_SHEET" > "${SAMPLE_SHEET}.tmp" && mv "${SAMPLE_SHEET}.tmp" "$SAMPLE_SHEET"

echo "Done! Removed duplicates from $SAMPLE_SHEET"
echo "Backup saved as: $BACKUP"
echo ""
echo "Verifying (last 12 lines):"
tail -12 "$SAMPLE_SHEET"


