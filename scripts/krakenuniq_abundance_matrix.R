#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  stop("Usage: Rscript krakenuniq_abundance_matrix.R <input_dir> <output_dir> <min_reads> <top_n>")
}
input_dir <- args[1]
output_dir <- args[2]
min_reads <- as.numeric(args[3])
top_n <- as.numeric(args[4])

# Find KrakenUniq report files in input_dir (look recursively)
files <- list.files(input_dir, pattern = "kraken-report.txt$", recursive = TRUE, full.names = TRUE)
cat("KrakenUniq report files found:\n")
print(files)

if(length(files) == 0) {
  stop("Error: No krakenuniq report files found in ", input_dir)
}

# Function to read one KrakenUniq report file and filter species >= min_reads
read_krakenuniq_report <- function(file, min_reads) {
  cat("Reading file:", file, "\n")

  empty_result <- tibble(
    reads = numeric(),
    rank = character(),
    taxID = character(),
    taxName = character(),
    sample = character()
  )
  
  if (!file.exists(file)) {
    warning("Skipping missing file: ", file, call. = FALSE)
    return(empty_result)
  }
  
  if (isTRUE(file.size(file) == 0)) {
    warning("Skipping empty KrakenUniq report: ", file, call. = FALSE)
    return(empty_result)
  }
  
  # Find the header line (starts with %) - robust to variable number of comment lines
  lines <- readLines(file, warn = FALSE)
  if (length(lines) == 0) {
    warning("Skipping empty KrakenUniq report (no lines): ", file, call. = FALSE)
    return(empty_result)
  }
  header_idx <- NULL
  for(i in seq_along(lines)) {
    # Check if line starts with % (possibly with leading whitespace)
    line_trimmed <- trimws(lines[i])
    if(startsWith(line_trimmed, "%") && grepl("\t", lines[i])) {
      header_idx <- i
      break
    }
  }
  
  if(is.null(header_idx)) {
    warning("Skipping KrakenUniq report with no header line: ", file, call. = FALSE)
    return(empty_result)
  }
  
  # Read the file starting from the header line
  # skip = header_idx - 1 means skip lines before header, then use header line as column names
  df <- readr::read_tsv(file, skip = header_idx - 1, show_col_types = FALSE, progress = FALSE)
  if (nrow(df) == 0) {
    return(empty_result)
  }

  # Normalize core column types across all input files so bind_rows() never fails
  # (readr may infer different types depending on file contents)
  df <- df %>%
    mutate(
      reads = suppressWarnings(as.numeric(reads)),
      rank = as.character(rank),
      taxID = as.character(taxID),
      taxName = as.character(taxName)
    )
  
  # Check if required columns exist
  required_cols <- c("reads", "rank", "taxID", "taxName")
  missing_cols <- setdiff(required_cols, colnames(df))
  if(length(missing_cols) > 0) {
    cat("Warning: Missing columns in", file, ":", paste(missing_cols, collapse = ", "), "\n")
    cat("Available columns:", paste(colnames(df), collapse = ", "), "\n")
    warning("Skipping KrakenUniq report with missing required columns: ", file, call. = FALSE)
    return(empty_result)
  }
  
  # Rename % column if it exists (it's optional, we don't actually use it)
  if("%" %in% colnames(df)) {
    df <- df %>% rename(percent_reads = `%`)
  }
  
  # Filter species-level taxa with reads >= min_reads
  df <- df %>%
    filter(rank == "species", reads >= min_reads)
  
  # Sample ID is in the report filename: {sample}_kraken-report.txt (layout-agnostic).
  sample_name <- sub("_kraken-report\\.txt$", "", basename(file), ignore.case = TRUE)
  if (!nzchar(sample_name)) {
    warning("Could not parse sample name from KrakenUniq report path: ", file, call. = FALSE)
    return(empty_result)
  }
  df <- df %>% mutate(sample = sample_name)
  return(df)
}

# Read all files and combine
all_species <- map_dfr(
  files,
  ~tryCatch(
    read_krakenuniq_report(.x, min_reads),
    error = function(e) {
      warning("Skipping KrakenUniq report due to read error: ", .x, "\n  ", conditionMessage(e), call. = FALSE)
      tibble(reads = numeric(), rank = character(), taxID = character(), taxName = character(), sample = character())
    }
  )
)

if(nrow(all_species) == 0) {
  stop("Error: No species-level taxa passed filtering. Check min_reads or input files.")
}

cat("Number of species-level taxa after filtering:", nrow(all_species), "\n")

# Summarize total reads per species across samples
species_totals <- all_species %>%
  group_by(taxID, taxName) %>%
  summarize(total_reads = sum(reads), .groups = "drop") %>%
  arrange(desc(total_reads))

# Select top N species by total_reads
top_species <- species_totals %>%
  slice_head(n = top_n) %>%
  pull(taxID)

# Filter to only top species
filtered_species <- all_species %>%
  filter(taxID %in% top_species)

# Create abundance matrix: rows = species, columns = samples
abundance_mat <- filtered_species %>%
  select(sample, taxName, reads) %>%
  pivot_wider(names_from = sample, values_from = reads, values_fill = 0)

# Set rownames to taxName and remove taxName column for matrix export
abundance_mat_df <- abundance_mat %>% column_to_rownames("taxName")

# Save absolute abundance matrix CSV
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(abundance_mat_df, file = file.path(output_dir, "krakenuniq_abundance_matrix_absolute.csv"), quote = FALSE)

# Create normalized abundance matrix (relative abundance per sample)
abundance_mat_norm <- sweep(as.matrix(abundance_mat_df), 2, colSums(as.matrix(abundance_mat_df)), FUN = "/")
abundance_mat_norm[is.na(abundance_mat_norm)] <- 0  # replace NaNs if any column sums zero

write.csv(abundance_mat_norm, file = file.path(output_dir, "krakenuniq_abundance_matrix_normalized.csv"), quote = FALSE)

cat("Abundance matrices saved in", output_dir, "\n")
