#!/usr/bin/env bash
# prepare_ena_dataset.sh
# Downloads ENA FASTQs from multiple filereports, creates sample sheet, and merges with local samples.
# Usage: bash prepare_ena_dataset.sh [--threads N] [--connections-per-file N] [--ena-dir DIR] [--ena-files FILE1 FILE2 ...]

set -eo pipefail

# ============================================================================
# Configuration (edit these as needed)
# ============================================================================

# Local samples file (if it exists, will be merged with ENA samples)
LOCAL_SAMPLES_TSV="${LOCAL_SAMPLES_TSV:-config/local_samples.tsv}"

# ENA filereport TSV files (glob pattern) - used if --ena-files not specified
ENA_FILEREPORT_GLOB="${ENA_FILEREPORT_GLOB:-config/ENA/filereport_read_run_*.tsv}"

# Where to download ENA FASTQs
ENA_OUTPUT_DIR="${ENA_OUTPUT_DIR:-/raid_md0/louis/ENA_raw_pigsti}"

# Output files
ENA_SAMPLES_TSV="config/ena_samples.tsv"
MERGED_SAMPLES_TSV="config/samples.tsv"

# Default number of parallel download threads
DEFAULT_THREADS=32

# Default aria2c connections per file
DEFAULT_CONNECTIONS_PER_FILE=2

# ============================================================================
# Parse command-line arguments
# ============================================================================

THREADS="${DEFAULT_THREADS}"
CONNECTIONS_PER_FILE="${DEFAULT_CONNECTIONS_PER_FILE}"
ENA_FILES_SPECIFIED=0
ENA_FILES_ARRAY=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --connections-per-file)
            CONNECTIONS_PER_FILE="$2"
            shift 2
            ;;
        --ena-dir)
            ENA_OUTPUT_DIR="$2"
            shift 2
            ;;
        --ena-files)
            # Collect all files until next option or end
            shift
            while [[ $# -gt 0 ]] && [[ ! "$1" =~ ^-- ]]; do
                ENA_FILES_ARRAY+=("$1")
                shift
            done
            ENA_FILES_SPECIFIED=1
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Downloads ENA FASTQs from multiple filereports and prepares merged samples.tsv for PIGSTI."
            echo ""
            echo "Options:"
            echo "  --threads N              Number of parallel download threads (default: ${DEFAULT_THREADS})"
            echo "  --connections-per-file N Number of aria2c connections per file (default: ${DEFAULT_CONNECTIONS_PER_FILE})"
            echo "  --ena-dir DIR            Directory for ENA FASTQ downloads (default: ${ENA_OUTPUT_DIR})"
            echo "  --ena-files FILE1 ...    Explicitly specify ENA filereport TSV files"
            echo "                          (if not specified, uses glob: ${ENA_FILEREPORT_GLOB})"
            echo "  --help, -h               Show this help message"
            echo ""
            echo "Examples:"
            echo "  # Use all files matching glob pattern"
            echo "  $0 --threads 32"
            echo ""
            echo "  # Specify files explicitly"
            echo "  $0 --threads 32 --ena-files config/ENA/filereport1.tsv config/ENA/filereport2.tsv"
            echo ""
            echo "Environment variables (optional):"
            echo "  LOCAL_SAMPLES_TSV      Path to local samples TSV (default: ${LOCAL_SAMPLES_TSV})"
            echo "  ENA_FILEREPORT_GLOB    Glob pattern for ENA filereports (default: ${ENA_FILEREPORT_GLOB})"
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            echo "Run '$0 --help' for usage information." >&2
            exit 1
            ;;
    esac
done

# ============================================================================
# Collect ENA filereport files
# ============================================================================

echo "============================================================================"
echo "PIGSTI ENA Dataset Preparation"
echo "============================================================================"
echo ""

if [ ${ENA_FILES_SPECIFIED} -eq 1 ]; then
    # Use explicitly specified files
    ENA_FILES=("${ENA_FILES_ARRAY[@]}")
    echo "Using explicitly specified ENA filereport files:"
else
    # Use glob pattern
    ENA_FILES=($(ls ${ENA_FILEREPORT_GLOB} 2>/dev/null || true))
    echo "Searching for ENA filereport files matching: ${ENA_FILEREPORT_GLOB}"
fi

# Validate that files exist
if [ ${#ENA_FILES[@]} -eq 0 ]; then
    echo "" >&2
    echo "ERROR: No ENA filereport files found!" >&2
    if [ ${ENA_FILES_SPECIFIED} -eq 1 ]; then
        echo "  Specified files:" >&2
        for f in "${ENA_FILES_ARRAY[@]}"; do
            echo "    - $f" >&2
        done
    else
        echo "  Glob pattern: ${ENA_FILEREPORT_GLOB}" >&2
        echo "  Please ensure ENA filereport TSV files exist in config/ENA/ directory." >&2
    fi
    exit 1
fi

# Check each file exists
MISSING_FILES=()
for f in "${ENA_FILES[@]}"; do
    if [ ! -f "$f" ]; then
        MISSING_FILES+=("$f")
    fi
done

if [ ${#MISSING_FILES[@]} -gt 0 ]; then
    echo "" >&2
    echo "ERROR: Some specified ENA filereport files do not exist:" >&2
    for f in "${MISSING_FILES[@]}"; do
        echo "  - $f" >&2
    done
    exit 1
fi

echo ""
echo "Found ${#ENA_FILES[@]} ENA filereport file(s):"
for f in "${ENA_FILES[@]}"; do
    echo "  - $f"
done
echo ""

# Check if Python script exists
if [ ! -f "scripts/ena_to_pigsti_samplesheet.py" ]; then
    echo "ERROR: scripts/ena_to_pigsti_samplesheet.py not found!" >&2
    echo "Please run this script from the PIGSTI root directory." >&2
    exit 1
fi

# ============================================================================
# Step 1: Download ENA FASTQs and create ENA sample sheet
# ============================================================================

echo "============================================================================"
echo "Step 1/3: Downloading ENA FASTQs (using ${THREADS} parallel workers, ${CONNECTIONS_PER_FILE} connections per file)"
echo "============================================================================"
echo ""

mkdir -p "$(dirname "${ENA_SAMPLES_TSV}")"
mkdir -p "${ENA_OUTPUT_DIR}"

python scripts/ena_to_pigsti_samplesheet.py \
    "${ENA_FILES[@]}" \
    --output-dir "${ENA_OUTPUT_DIR}" \
    --samples-tsv "${ENA_SAMPLES_TSV}" \
    --max-workers "${THREADS}" \
    --connections-per-file "${CONNECTIONS_PER_FILE}"

if [ ! -f "${ENA_SAMPLES_TSV}" ]; then
    echo "ERROR: ENA sample sheet was not created: ${ENA_SAMPLES_TSV}" >&2
    exit 1
fi

echo ""
echo "✓ ENA downloads complete. Sample sheet: ${ENA_SAMPLES_TSV}"
echo ""

# ============================================================================
# Step 2: Prepare local samples (if they exist)
# ============================================================================

echo "============================================================================"
echo "Step 2/3: Preparing local samples"
echo "============================================================================"
echo ""

LOCAL_TEMP=""
BASE_LOCAL_TSV=""

# Prefer explicit local_samples.tsv; fall back to an existing merged samples.tsv
if [ -f "${LOCAL_SAMPLES_TSV}" ]; then
    BASE_LOCAL_TSV="${LOCAL_SAMPLES_TSV}"
    echo "Found local samples file: ${BASE_LOCAL_TSV}"
elif [ -f "${MERGED_SAMPLES_TSV}" ]; then
    BASE_LOCAL_TSV="${MERGED_SAMPLES_TSV}"
    echo "Found existing samples sheet, treating as local: ${BASE_LOCAL_TSV}"
fi

if [ -n "${BASE_LOCAL_TSV}" ]; then
    # Check if local samples already have 'source' column
    if head -n 1 "${BASE_LOCAL_TSV}" | grep -q "source"; then
        echo "  Local samples already have 'source' column, using as-is."
        LOCAL_TEMP="${BASE_LOCAL_TSV}"
    else
        echo "  Adding 'source=LOCAL' column to local samples."
        LOCAL_TEMP=$(mktemp)
        # Add 'source' column header and set all rows to 'LOCAL'
        head -n 1 "${BASE_LOCAL_TSV}" | awk -F'\t' '{print $0 "\tsource"}' > "${LOCAL_TEMP}"
        tail -n +2 "${BASE_LOCAL_TSV}" | awk -F'\t' '{print $0 "\tLOCAL"}' >> "${LOCAL_TEMP}"
    fi
else
    echo "No local samples file found: ${LOCAL_SAMPLES_TSV}"
    echo "No existing samples sheet found at ${MERGED_SAMPLES_TSV}"
    echo "  (Skipping local samples merge)"
fi
echo ""

# ============================================================================
# Step 3: Merge local + ENA into final samples.tsv
# ============================================================================

echo "============================================================================"
echo "Step 3/3: Merging local + ENA samples into ${MERGED_SAMPLES_TSV}"
echo "============================================================================"
echo ""

mkdir -p "$(dirname "${MERGED_SAMPLES_TSV}")"

# Get header from ENA samples (it should have 'source' column)
HEADER=$(head -n 1 "${ENA_SAMPLES_TSV}")

# Start merged file with header
echo "${HEADER}" > "${MERGED_SAMPLES_TSV}"

# Append local samples (skip header if using temp file, or skip header from original)
if [ -n "${LOCAL_TEMP}" ]; then
    tail -n +2 "${LOCAL_TEMP}" >> "${MERGED_SAMPLES_TSV}"
    echo "  Added $(tail -n +2 "${LOCAL_TEMP}" | wc -l) local sample(s)"
fi

# Append ENA samples (skip header)
tail -n +2 "${ENA_SAMPLES_TSV}" >> "${MERGED_SAMPLES_TSV}"
ENA_COUNT=$(tail -n +2 "${ENA_SAMPLES_TSV}" | wc -l)
echo "  Added ${ENA_COUNT} ENA sample(s)"

# Clean up temp file if we created one
if [ -n "${LOCAL_TEMP}" ] && [ "${LOCAL_TEMP}" != "${LOCAL_SAMPLES_TSV}" ]; then
    rm -f "${LOCAL_TEMP}"
fi

TOTAL_COUNT=$(tail -n +2 "${MERGED_SAMPLES_TSV}" | wc -l)
echo ""
echo "✓ Merged sample sheet created: ${MERGED_SAMPLES_TSV}"
echo "  Total samples: ${TOTAL_COUNT}"
echo ""

# ============================================================================
# Summary
# ============================================================================

echo "============================================================================"
echo "✓ ENA Dataset Preparation Complete!"
echo "============================================================================"
echo ""
echo "Next steps:"
echo "  1. Review the merged sample sheet: ${MERGED_SAMPLES_TSV}"
echo "  2. Run PIGSTI pipeline:"
echo "     snakemake --use-conda --cores 42"
echo ""
echo "Note: Raw FASTQs are kept on disk (PIGSTI does not delete input FASTQ files)."
echo "      Raw read counts are stored under results/libraries/{pcr}/qc/*.raw_reads.txt"
echo ""

