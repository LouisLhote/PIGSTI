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
  df <- readr::read_tsv(file, skip = 2, show_col_types = FALSE) %>%
    rename(percent_reads = `%`) %>%
    filter(rank == "species", reads >= min_reads)
  # Add sample name from file path (take first folder under input_dir)
  sample_name <- str_split(file, .Platform$file.sep)[[1]][length(str_split(input_dir, .Platform$file.sep)[[1]]) + 1]
  df <- df %>% mutate(sample = sample_name)
  return(df)
}

# Read all files and combine
all_species <- map_dfr(files, ~read_krakenuniq_report(.x, min_reads))

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
