#!/usr/bin/env Rscript

# ------------------------- LIBRARIES -------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(patchwork)
  library(pdftools)
  library(magick)
  library(readr)
  library(viridis)
  library(grid)
})

# ------------------------- INPUT ARGS -------------------------
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
pathogen <- args[2]
spreadsheet_path <- args[3]
output_path <- args[4]

safe_pathogen <- gsub(" ", "_", pathogen)
results_dir <- file.path("results", sample)

# ------------------------- FILE PATHS -------------------------
escore_file <- file.path(results_dir, "Escore/pathogen", paste0(sample, "_pathogen.csv"))
ani_file <- file.path(results_dir, "bwa_pathogen", paste0(sample, "_", safe_pathogen, ".ani.txt"))
entropy_plot <- file.path(results_dir, "bwa_pathogen", paste0(sample, "_", safe_pathogen, ".entropy_plot.png"))
bamplot_pdf <- file.path(results_dir, "bwa_pathogen", paste0("adnaplotter_", safe_pathogen, ".pdf"))
damage_dir <- file.path(results_dir, "bwa_pathogen", paste0("damageprofiler_", safe_pathogen))
damage_pdfs <- list(
  mis = file.path(damage_dir, "DamagePlot.pdf"),
  edit = file.path(damage_dir, "edit_distance.pdf"),
  length = file.path(damage_dir, "Length_plot.pdf")
)

# ------------------------- SPREADSHEET -------------------------
spreadsheet <- read.csv(spreadsheet_path)
spreadsheet$Krakenuniq.name <- trimws(spreadsheet$Krakenuniq.name)
spreadsheet$Hops.name <- trimws(spreadsheet$Hops.name)

# ------------------------- HOPS SUMMARY PDF -------------------------
# ------------------------- HOPS REPORT PDF -------------------------
summary_pdf <- NULL
hops_name <- spreadsheet$Hops.name[spreadsheet$Krakenuniq.name == pathogen]

if (length(hops_name) > 0 && !is.na(hops_name)) {
  hops_name <- trimws(hops_name)
  safe_hops_name <- gsub(" ", "_", hops_name)

  hops_dir <- file.path("results", "hops", "maltExtract", "pdf_candidate_profiles", safe_hops_name)

  # Use pattern matching to locate file: stpX_SAMPLE_unaligned.rma6_HOPSNAME_summary.pdf
  summary_files <- list.files(
    path = hops_dir,
    pattern = paste0("^stp[0-9]+_", sample, "_unaligned\\.rma6_", safe_hops_name, "_summary\\.pdf$"),
    full.names = TRUE
  )

  if (length(summary_files) > 0 && file.exists(summary_files[1])) {
    summary_pdf <- summary_files[1]
  } else {
    message("[WARNING] No HOPS summary PDF found in: ", hops_dir)
  }
}

# ------------------------- PLOT FUNCTIONS -------------------------
plot_escore <- function() {
  if (!file.exists(escore_file)) return(NULL)
  escore <- read.csv(escore_file)
  escore$taxonomy <- trimws(escore$taxonomy)
  escore <- escore[escore$taxonomy == pathogen, ]

  if (nrow(escore) == 0) return(NULL)
  row <- spreadsheet[spreadsheet$Krakenuniq.name == pathogen, ]
  threshold_escore <- row$min_escore
  threshold_reads <- row$min_reads

  ggplot(escore, aes(x = Escore, y = taxReads)) +
    geom_point(color = "steelblue", size = 3) +
    geom_vline(xintercept = threshold_escore, linetype = "dashed", color = "red") +
    geom_hline(yintercept = threshold_reads, linetype = "dashed", color = "red") +
    labs(title = paste("Escore vs Reads for", pathogen),
         x = "Escore", y = "Read Count") +
    theme_minimal()
}

plot_ani <- function() {
  if (!file.exists(ani_file)) return(NULL)
  ani_val <- as.numeric(gsub("[^0-9\\.]", "", readLines(ani_file)))
  if (is.na(ani_val)) return(NULL)
  df <- data.frame(ANI = ani_val)
  ggplot(df, aes(x = "", y = ANI)) +
    geom_col(fill = "forestgreen", width = 0.5) +
    coord_flip() +
    ylim(0, 100) +
    labs(title = "ANI", y = "% Identity", x = NULL) +
    theme_minimal()
}

plot_entropy <- function() {
  if (!file.exists(entropy_plot)) return(NULL)
  img <- image_read(entropy_plot)
  rasterGrob(img, interpolate = TRUE)
}

embed_pdf_page <- function(pdf_path) {
  if (!file.exists(pdf_path)) return(NULL)
  img <- image_read_pdf(pdf_path, density = 150)
  rasterGrob(img, interpolate = TRUE)
}

# ------------------------- LAYOUT AND EXPORT -------------------------
panels <- list(
  A1 = plot_escore(),
  A2 = if (!is.null(summary_pdf) && file.exists(summary_pdf)) embed_pdf_page(summary_pdf) else NULL,
  A3 = NULL,
  B1 = plot_entropy(),
  B2 = plot_ani(),
  B3 = embed_pdf_page(bamplot_pdf),
  C1 = embed_pdf_page(damage_pdfs$mis),
  C2 = embed_pdf_page(damage_pdfs$edit),
  C3 = embed_pdf_page(damage_pdfs$length)
)

valid_panels <- panels[!sapply(panels, is.null)]

if (length(valid_panels) == 0) {
  message("[WARNING] No plots available for ", sample, " - ", pathogen)
  quit(status = 0)
}

layout <- marrangeGrob(grobs = valid_panels, ncol = 3, nrow = 3, top = paste(sample, "-", pathogen))
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
ggsave(output_path, layout, width = 11.7, height = 8.3)  # A4 landscape
