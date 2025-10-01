#!/usr/bin/env Rscript

if (!requireNamespace("MetBrewer", quietly = TRUE)) {
  install.packages("MetBrewer", repos = "https://cloud.r-project.org")
}

library(tidyverse)
library(MetBrewer)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript krakenuniq_plot_heatmap.R <indir> <abs_plot_out.pdf> <norm_plot_out.pdf>")
}

indir <- args[1]               # e.g. "results/KRAKENUNIQ_ABUNDANCE_MATRIX"
abs_plot_out <- args[2]        # PDF output for absolute heatmap
norm_plot_out <- args[3]       # PDF output for normalized heatmap

# Load matrices
abs_path <- file.path(indir, "krakenuniq_abundance_matrix_absolute.csv")
norm_path <- file.path(indir, "krakenuniq_abundance_matrix_normalized.csv")

df_abs <- read.csv(abs_path, row.names = 1, check.names = FALSE)
df_norm <- read.csv(norm_path, row.names = 1, check.names = FALSE)

# Convert to long format
df_abs_long <- df_abs %>%
  rownames_to_column("species") %>%
  pivot_longer(-species, names_to = "sample", values_to = "abundance")

df_norm_long <- df_norm %>%
  rownames_to_column("species") %>%
  pivot_longer(-species, names_to = "sample", values_to = "norm_abundance") %>%
  mutate(norm_label = scales::percent(norm_abundance, accuracy = 0.1))

# Order species by total absolute abundance
species_order <- df_abs_long %>%
  group_by(species) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  pull(species)

df_abs_long <- df_abs_long %>%
  mutate(species = factor(species, levels = species_order))

df_norm_long <- df_norm_long %>%
  mutate(species = factor(species, levels = species_order))

# Define color palette
hokusai_pal <- met.brewer("Hokusai1", type = "continuous", n = 50)

# -------- Absolute Abundance Heatmap --------
abs_plot <- ggplot(df_abs_long, aes(x = sample, y = species, fill = abundance)) +
  geom_tile() +
  geom_text(aes(label = abundance), size = 3, color = "white") +
  scale_fill_gradientn(colors = hokusai_pal, na.value = "grey90") +
  theme_minimal() +
  labs(title = "Absolute Abundance Heatmap", x = "Sample", y = "Species") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold")
  )

ggsave(abs_plot_out, abs_plot, width = 14, height = 10)

# -------- Normalized Abundance Heatmap --------
norm_plot <- ggplot(df_norm_long, aes(x = sample, y = species, fill = norm_abundance)) +
  geom_tile() +
  geom_text(aes(label = norm_label), size = 3, color = "white") +
  scale_fill_gradientn(colors = hokusai_pal, na.value = "grey90") +
  theme_minimal() +
  labs(title = "Normalized Abundance Heatmap", x = "Sample", y = "Species") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold")
  )

ggsave(norm_plot_out, norm_plot, width = 14, height = 10)

cat("Heatmaps saved to:\n")
cat("  Absolute:", abs_plot_out, "\n")
cat("  Normalized:", norm_plot_out, "\n")
