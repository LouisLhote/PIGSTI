#!/usr/bin/env Rscript

# ------------------------- LIBRARIES -------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(patchwork)
  library(pdftools)
  library(magick)
  library(readr)
  library(viridis)
  library(ggthemes)
  library(MetBrewer)
  library(scales)
  library(jsonlite)
})

# Try to load extrafont for Lato font
tryCatch({
  if (requireNamespace("extrafont", quietly = TRUE)) {
    library(extrafont)
    extrafont::loadfonts(device = "pdf", quiet = TRUE)
  }
}, error = function(e) {
  message("[INFO] extrafont not available, using default fonts")
})

# ------------------------- INPUT ARGS -------------------------
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
pathogen <- args[2]
spreadsheet_path <- args[3]
output_path <- args[4]

safe_pathogen <- gsub(" ", "_", pathogen)
results_dir <- file.path("results", "pathogen", sample)

# ------------------------- FILE PATHS -------------------------
escore_file <- file.path(results_dir, "evalue", "pathogen", paste0(sample, "_pathogen.csv"))
summary_file <- file.path(results_dir, "summary", paste0(sample, "_pathogen_summary.csv"))
ani_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".ani.txt"))
ani_dist_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".ani_distribution.txt"))
entropy_plot <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".entropy_plot.png"))
entropy_mean_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".mean_entropy.txt"))
# AdnaPlotter removed - too time and memory consuming
damage_dir <- file.path(results_dir, "pathogen_mapping", paste0("damageprofiler_", safe_pathogen))
qualimap_file <- file.path(results_dir, "pathogen_mapping", paste0("qualimap_", safe_pathogen), "genome_results.txt")
breadth_ratio_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".breadth_ratio.txt"))
depth_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".depth.txt"))
bam_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".dedup.bam"))
edit_r2_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".edit_distance_logr2.txt"))
edit_r2_plot <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".edit_distance_logr2.png"))
edit_r2_dist_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".edit_distance_logr2_distribution.txt"))

# Optional HOPS-like damage-vs-no-damage edit distance metrics
edit_r2_damaged_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".edit_distance_logr2_damaged.txt"))
edit_r2_damaged_dist_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".edit_distance_logr2_damaged_distribution.txt"))
edit_r2_no_damage_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".edit_distance_logr2_no_damage.txt"))
edit_r2_no_damage_dist_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".edit_distance_logr2_no_damage_distribution.txt"))
genus_ranking_file <- file.path(results_dir, "pathogen_mapping", paste0(sample, "_", safe_pathogen, ".genus_ranking.txt"))
kraken_report_file <- file.path(results_dir, "krakenuniq", paste0(sample, "_kraken-report.txt"))

damage_pdfs <- list(
  mis = file.path(damage_dir, "DamagePlot.pdf"),
  edit = file.path(damage_dir, "edit_distance.pdf"),
  length = file.path(damage_dir, "Length_plot.pdf")
)

# ------------------------- LOAD DETECTION CRITERIA CONFIG -------------------------
criteria_config_path <- file.path(dirname(spreadsheet_path), "pathogen_detection_criteria.yaml")
criteria_config <- list(
  ani_threshold = 96.5,
  entropy_threshold = 0.9,
  entropy_threshold_virus = 0.7,
  breadth_ratio_threshold = 0.8,
  damage_5p_ct_threshold = 0.01,
  damage_3p_ga_threshold = 0.01,
  hops_detection_threshold = 2,
  hops_edit_distance_threshold = 3,
  hops_damage_threshold = 4,
  high_confidence_threshold = 0.7,
  edit_distance_decay_quality_threshold = 0.65,
  edit_distance_logr2_threshold = 0.60,
  read_mapping_ratio_threshold = 0.5,
  escore_threshold = 5,
  reads_threshold = 50,
  guellil_evalue_threshold = 0.001,
  use_evalue_for_detection = TRUE,
  edit_distance_no_damage_decay_quality_min = 0.55,
  edit_distance_damage_minus_no_damage_min = 0.10
)

# Prefer embedded criteria in config/config.yaml (pathogen_detection_criteria)
main_config_path <- "config/config.yaml"
if (file.exists(main_config_path) && requireNamespace("yaml", quietly = TRUE)) {
  tryCatch({
    main_cfg <- yaml::read_yaml(main_config_path)
    if (!is.null(main_cfg$pathogen_detection_criteria)) {
      criteria_config <- modifyList(criteria_config, main_cfg$pathogen_detection_criteria)
    }
  }, error = function(e) {
    message("[WARNING] Could not read config/config.yaml pathogen_detection_criteria: ", e$message)
  })
}

if (file.exists(criteria_config_path)) {
  tryCatch({
    if (requireNamespace("yaml", quietly = TRUE)) {
      loaded_cfg <- yaml::read_yaml(criteria_config_path)
      # Merge loaded config on top of defaults, so missing keys keep defaults.
      criteria_config <- modifyList(criteria_config, loaded_cfg)
    } else {
      # Fallback: try to parse as simple key-value pairs
      config_lines <- readLines(criteria_config_path)
      for (line in config_lines) {
        if (grepl("^[^#].*:", line)) {
          parts <- strsplit(line, ":")[[1]]
          if (length(parts) == 2) {
            key <- trimws(parts[1])
            value <- trimws(parts[2])
            # Try to convert to numeric
            num_value <- suppressWarnings(as.numeric(value))
            if (!is.na(num_value)) {
              criteria_config[[key]] <- num_value
            }
          }
        }
      }
    }
  }, error = function(e) {
    message("[WARNING] Could not load criteria config, using defaults: ", e$message)
  })
}

# ------------------------- HELPER FUNCTIONS -------------------------
safe_read_csv <- function(path, ..., allow_empty = TRUE) {
  if (is.null(path) || is.na(path) || path == "") return(NULL)
  if (!file.exists(path)) return(NULL)
  if (isTRUE(file.size(path) == 0)) {
    if (!allow_empty) stop("CSV file is empty: ", path)
    return(NULL)
  }
  first_line <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) character(0))
  if (length(first_line) == 0) {
    if (!allow_empty) stop("CSV file has no readable lines: ", path)
    return(NULL)
  }
  tryCatch(
    read.csv(path, ...),
    error = function(e) {
      if (!allow_empty) stop("Failed to read CSV: ", path, "\n", conditionMessage(e))
      message("[WARNING] Skipping unreadable CSV: ", path, " (", conditionMessage(e), ")")
      NULL
    }
  )
}

# ------------------------- SPREADSHEET -------------------------
spreadsheet <- safe_read_csv(spreadsheet_path, allow_empty = FALSE)
spreadsheet$Krakenuniq.name <- trimws(spreadsheet$Krakenuniq.name)
spreadsheet$Hops.name <- trimws(spreadsheet$Hops.name)

# ------------------------- READ SUMMARY DATA -------------------------
summary_data <- NULL
summary_data_raw <- safe_read_csv(summary_file, allow_empty = TRUE)
if (!is.null(summary_data_raw) && nrow(summary_data_raw) > 0) {
  if ("Pathogen" %in% colnames(summary_data_raw)) {
    summary_data <- summary_data_raw[summary_data_raw$Pathogen == pathogen, ]
  } else {
    message("[WARNING] Summary CSV missing 'Pathogen' column: ", summary_file)
    summary_data <- NULL
  }
}

# ------------------------- HOPS SUMMARY PDF -------------------------
summary_pdf <- NULL
hops_name <- spreadsheet$Hops.name[spreadsheet$Krakenuniq.name == pathogen]

if (length(hops_name) > 0 && !is.na(hops_name)) {
  hops_name <- trimws(hops_name)
  safe_hops_name <- gsub(" ", "_", hops_name)

  hops_dir <- file.path("results", "hops", "maltExtract", "pdf_candidate_profiles", safe_hops_name)

  if (dir.exists(hops_dir)) {
  summary_files <- list.files(
    path = hops_dir,
    pattern = paste0("^stp[0-9]+_", sample, "_unaligned\\.rma6_", safe_hops_name, "_summary\\.pdf$"),
    full.names = TRUE
  )

  if (length(summary_files) > 0 && file.exists(summary_files[1])) {
    summary_pdf <- summary_files[1]
    }
  }
}

read_metric_file <- function(file_path) {
  if (!file.exists(file_path)) return(NA)
  tryCatch({
    val <- as.numeric(readLines(file_path, n = 1))
    return(val)
  }, error = function(e) NA)
}

extract_qualimap_coverage <- function(file_path) {
  if (!file.exists(file_path)) return(NA)
  tryCatch({
    lines <- readLines(file_path)
    for (line in lines) {
      if (grepl("mean coverage", line, ignore.case = TRUE)) {
        val_str <- strsplit(line, "=")[[1]][2]
        val <- suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", val_str)))
        if (!is.na(val)) return(val)
      }
    }
    return(NA)
  }, error = function(e) NA)
}

extract_qualimap_evenness <- function(file_path) {
  if (!file.exists(file_path)) return(NA)
  tryCatch({
    lines <- readLines(file_path)
    for (line in lines) {
      # Qualimap line looks like:
      # "There is a 12.345% of reference with a coverageData >= 1X"
      if (grepl("coverageData\\s*>=\\s*1X", line, ignore.case = TRUE)) {
        percent_str <- strsplit(line, "There is a", fixed = TRUE)[[1]]
        if (length(percent_str) > 1) {
          pct <- trimws(strsplit(percent_str[2], "%")[[1]][1])
          val <- suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", pct)))
          if (!is.na(val)) return(val)
        }
      }
    }
    return(NA)
  }, error = function(e) NA)
}

read_escore_row <- function(escore_path, pathogen_name) {
  if (!file.exists(escore_path)) return(NULL)
  esc <- safe_read_csv(escore_path, allow_empty = TRUE)
  if (is.null(esc) || nrow(esc) == 0) return(NULL)
  if (!"taxonomy" %in% colnames(esc)) return(NULL)
  esc$taxonomy <- trimws(as.character(esc$taxonomy))
  hit <- esc[esc$taxonomy == pathogen_name, , drop = FALSE]
  if (nrow(hit) == 0) return(NULL)
  return(hit[1, , drop = FALSE])
}

extract_mapped_reads_from_bam <- function(bam_path) {
  if (!file.exists(bam_path)) return(NA)
  # Use Rsamtools if available (fast, uses BAM index where possible)
  if (requireNamespace("Rsamtools", quietly = TRUE)) {
    tryCatch({
      stats <- Rsamtools::idxstatsBam(bam_path)
      mapped <- sum(stats$mapped, na.rm = TRUE)
      if (is.finite(mapped)) return(mapped)
      return(NA)
    }, error = function(e) NA)
  }
  return(NA)
}

# ------------------------- BEAUTIFUL THEME -------------------------
# Try to load Lato font, fallback to sans if not available
font_family <- "sans"  # Default fallback

tryCatch({
  if (requireNamespace("extrafont", quietly = TRUE)) {
    # Try to load fonts
    extrafont::loadfonts(device = "pdf", quiet = TRUE)
    
    # Check if Lato is actually available
    available_fonts <- extrafont::fonts()
    if ("Lato" %in% available_fonts) {
      font_family <- "Lato"
    } else {
      font_family <- "sans"
      message("[INFO] Lato font not found, using default sans font")
    }
  } else {
    font_family <- "sans"
    message("[INFO] extrafont package not available, using default sans font")
  }
}, error = function(e) {
  font_family <<- "sans"
  message("[INFO] Error loading fonts, using default sans font: ", e$message)
})

beautiful_theme <- function(base_size = 16) {
  # Only specify font family if it's not the default "sans"
  if (font_family != "sans") {
    theme_few(base_size = base_size, base_family = font_family) +
      theme(
        plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5, 
                                  margin = margin(b = 15), family = font_family),
        plot.subtitle = element_text(size = base_size - 2, hjust = 0.5, color = "gray40", family = font_family),
        axis.title = element_text(size = base_size, face = "bold", family = font_family),
        axis.text = element_text(size = base_size - 2, family = font_family),
        legend.title = element_text(size = base_size, face = "bold", family = font_family),
        legend.text = element_text(size = base_size - 2, family = font_family),
        plot.margin = margin(20, 20, 20, 20),
        text = element_text(family = font_family)
      )
  } else {
    # Use default sans font (don't specify family to avoid PostScript errors)
    theme_few(base_size = base_size) +
      theme(
        plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5, 
                                  margin = margin(b = 15)),
        plot.subtitle = element_text(size = base_size - 2, hjust = 0.5, color = "gray40"),
        axis.title = element_text(size = base_size, face = "bold"),
        axis.text = element_text(size = base_size - 2),
        legend.title = element_text(size = base_size, face = "bold"),
        legend.text = element_text(size = base_size - 2),
        plot.margin = margin(20, 20, 20, 20)
      )
  }
}

# Beautiful color palettes from MetBrewer - use different palettes for different plots
get_met_palette <- function(palette_name = "VanGogh3", n = 8) {
  tryCatch({
    MetBrewer::met.brewer(palette_name, n = n)
  }, error = function(e) {
    # Fallback to beautiful colors if MetBrewer fails
    c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB", "#FF6B6B")
  })
}

# Default palette (for backward compatibility)
met_palette <- get_met_palette("VanGogh3", 8)

# ------------------------- TITLE PAGE -------------------------
create_title_page <- function() {
  # Extract metrics
  coverage <- if (!is.null(summary_data) && nrow(summary_data) > 0) {
    cov_str <- summary_data$Coverage[1]
    if (!is.na(cov_str) && cov_str != "NA") {
      as.numeric(gsub("[^0-9\\.]", "", cov_str))
    } else {
      extract_qualimap_coverage(qualimap_file)
    }
  } else {
    extract_qualimap_coverage(qualimap_file)
  }
  
  mapped_reads <- if (!is.null(summary_data) && nrow(summary_data) > 0) {
    mr <- summary_data$Mapped_reads[1]
    # Handle both numeric and character values like "20724" or "20724.0"
    if (is.numeric(mr)) {
      if (!is.na(mr)) mr else NA
    } else {
      suppressWarnings(as.numeric(as.character(mr)))
    }
  } else NA

  # Fallback: derive mapped reads from BAM if summary is missing/empty
  if (is.na(mapped_reads)) {
    mapped_reads <- extract_mapped_reads_from_bam(bam_file)
  }
  
  kraken_reads <- if (!is.null(summary_data) && nrow(summary_data) > 0) {
    kr <- summary_data$Krakenuniq_reads[1]
    if (is.numeric(kr)) {
      if (!is.na(kr)) kr else NA
    } else {
      suppressWarnings(as.numeric(as.character(kr)))
    }
  } else NA

  # Fallback: get KrakenUniq reads from Escore table if summary is missing/empty
  if (is.na(kraken_reads) && file.exists(escore_file)) {
    esc_row <- read_escore_row(escore_file, pathogen)
    if (!is.null(esc_row)) {
      # Prefer 'reads' (clade reads); fall back to 'taxReads' if needed.
      if ("reads" %in% colnames(esc_row)) {
        kraken_reads <- suppressWarnings(as.numeric(as.character(esc_row$reads[1])))
      } else if ("taxReads" %in% colnames(esc_row)) {
        kraken_reads <- suppressWarnings(as.numeric(as.character(esc_row$taxReads[1])))
      }
    }
  }
  
  # Prefer the pre-computed Read_mapping_ratio column if it exists and is numeric;
  # otherwise, compute it on the fly from Mapped_reads / Krakenuniq_reads
  read_mapping_ratio <- NA
  if (!is.null(summary_data) && nrow(summary_data) > 0 && "Read_mapping_ratio" %in% colnames(summary_data)) {
    rmr_raw <- summary_data$Read_mapping_ratio[1]
    if (is.numeric(rmr_raw)) {
      if (!is.na(rmr_raw)) {
        read_mapping_ratio <- rmr_raw
      }
    } else if (!is.na(rmr_raw) && rmr_raw != "") {
      rmr_num <- suppressWarnings(as.numeric(as.character(rmr_raw)))
      if (!is.na(rmr_num)) {
        read_mapping_ratio <- rmr_num
      }
    }
  }
  # If we still don't have a ratio, derive it from mapped_reads and kraken_reads
  if ((is.na(read_mapping_ratio) || !is.finite(read_mapping_ratio)) &&
      !is.na(mapped_reads) && mapped_reads > 0 &&
      !is.na(kraken_reads) && kraken_reads > 0) {
    read_mapping_ratio <- mapped_reads / kraken_reads
  }
  
  evenness <- if (!is.null(summary_data) && nrow(summary_data) > 0) {
    ev <- summary_data$Evenness[1]
    if (!is.na(ev) && ev != "NA") {
      as.numeric(gsub("[^0-9\\.]", "", ev))
    } else NA
  } else NA

  # Fallback: derive evenness from Qualimap genome_results.txt when available
  if (is.na(evenness)) {
    evenness <- extract_qualimap_evenness(qualimap_file)
  }

  # Prefer ANI from summary table (already parsed in summarize_pathogen_data),
  # fall back to raw .ani.txt metric file if needed
  ani <- if (!is.null(summary_data) && nrow(summary_data) > 0 && "ANI" %in% colnames(summary_data)) {
    ani_str <- summary_data$ANI[1]
    if (!is.na(ani_str) && ani_str != "NA") {
      suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", as.character(ani_str))))
    } else NA
  } else NA
  if (is.na(ani)) {
    ani_raw <- read_metric_file(ani_file)
    if (!is.na(ani_raw)) {
      ani <- suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", as.character(ani_raw))))
    }
  }
  
  breadth_ratio <- read_metric_file(breadth_ratio_file)
  
  entropy_mean <- read_metric_file(entropy_mean_file)
  
  damage_5p <- if (!is.null(summary_data) && nrow(summary_data) > 0) {
    dmg <- summary_data$Damage_5p_CtoT[1]
    if (is.numeric(dmg) && !is.na(dmg)) dmg else NA
  } else NA
  
  # Create beautiful title page
  title_grob <- grobTree(
    # Background with solid color (R doesn't support CSS gradients)
    rectGrob(gp = gpar(fill = "#f0f4f8", alpha = 0.3)),
    # Title
    textGrob("Pathogen Analysis Report", 
             gp = gpar(fontsize = 32, fontface = "bold", col = "#2c3e50"),
             y = 0.88, x = 0.5),
    # Sample and Pathogen
    textGrob(paste("Sample:", sample), 
             gp = gpar(fontsize = 20, fontface = "bold", col = "#34495e"),
             y = 0.75, x = 0.5),
    textGrob(paste("Pathogen:", pathogen), 
             gp = gpar(fontsize = 20, fontface = "bold", col = "#34495e"),
             y = 0.68, x = 0.5),
    # Section header
    textGrob("Summary Metrics", 
             gp = gpar(fontsize = 18, fontface = "bold", col = "#2c3e50"),
             y = 0.58, x = 0.5),
    # Metrics in two columns (larger font)
    textGrob(paste("Coverage:", ifelse(is.na(coverage), "N/A", paste0(round(coverage, 2), "X"))),
             gp = gpar(fontsize = 18, col = "#34495e"),
             y = 0.50, just = "left", x = 0.15),
    textGrob(paste("Mapped Reads:", ifelse(is.na(mapped_reads), "N/A", format(mapped_reads, big.mark = ","))),
             gp = gpar(fontsize = 18, col = "#34495e"),
             y = 0.50, just = "left", x = 0.60),
    textGrob(paste("KrakenUniq Reads:", ifelse(is.na(kraken_reads), "N/A", format(kraken_reads, big.mark = ","))),
             gp = gpar(fontsize = 18, col = "#34495e"),
             y = 0.43, just = "left", x = 0.15),
    textGrob(paste("Read Mapping Ratio:", ifelse(is.na(read_mapping_ratio), "N/A", 
             round(read_mapping_ratio, 3))),
             gp = gpar(fontsize = 18, col = "#34495e"),
             y = 0.43, just = "left", x = 0.60),
    textGrob(paste("ANI:", ifelse(is.na(ani), "N/A", paste0(round(ani, 2), "%"))),
             gp = gpar(fontsize = 18, col = "#34495e"),
             y = 0.36, just = "left", x = 0.15),
    textGrob(paste("Breadth Ratio:", ifelse(is.na(breadth_ratio), "N/A", paste0(round(breadth_ratio, 3)))),
             gp = gpar(fontsize = 18, col = "#34495e"),
             y = 0.36, just = "left", x = 0.60),
    textGrob(paste("Evenness:", ifelse(is.na(evenness), "N/A", paste0(round(evenness, 2), "%"))),
             gp = gpar(fontsize = 18, col = "#34495e"),
             y = 0.29, just = "left", x = 0.15),
    textGrob(paste("5' C>T Damage:", ifelse(is.na(damage_5p), "N/A", paste0(round(damage_5p, 4)))),
             gp = gpar(fontsize = 18, col = "#34495e"),
             y = 0.29, just = "left", x = 0.60),
    textGrob(paste("Mean Entropy:", ifelse(is.na(entropy_mean), "N/A", paste0(round(entropy_mean, 4)))),
             gp = gpar(fontsize = 18, col = "#34495e"),
             y = 0.22, just = "left", x = 0.15)
  )
  
  return(title_grob)
}

# ------------------------- PLOT FUNCTIONS -------------------------
plot_escore <- function() {
  if (!file.exists(escore_file)) return(NULL)
  escore <- safe_read_csv(escore_file, allow_empty = TRUE)
  if (is.null(escore) || nrow(escore) == 0) return(NULL)
  escore$taxonomy <- trimws(escore$taxonomy)
  escore <- escore[escore$taxonomy == pathogen, ]

  if (nrow(escore) == 0) return(NULL)

  # Decide which KrakenUniq read-count column to use, to match summarize_pathogen_data.py
  # KrakenUniq report header is typically:
  # "%, reads, taxReads, kmers, dup, cov, taxID, rank, taxName"
  #  - 'reads'    = reads in clade (including children)
  #  - 'taxReads' = reads assigned directly to this taxon
  read_col <- NULL
  if ("reads" %in% colnames(escore)) {
    read_col <- "reads"
  } else if ("taxReads" %in% colnames(escore)) {
    read_col <- "taxReads"
  } else if ("X..reads" %in% colnames(escore)) {
    # Sometimes '%' column name becomes 'X..' etc. after read.csv; handle common variants
    read_col <- "X..reads"
  } else if ("X.reads" %in% colnames(escore)) {
    read_col <- "X.reads"
  } else if ("X.of.reads" %in% colnames(escore)) {
    read_col <- "X.of.reads"
  } else {
    # Fallback: if no explicit reads column found, try the first numeric column after taxonomy
    numeric_cols <- setdiff(colnames(escore), c("taxonomy"))
    numeric_cols <- numeric_cols[sapply(escore[numeric_cols], is.numeric)]
    if (length(numeric_cols) > 0) {
      read_col <- numeric_cols[1]
    } else {
      warning("Could not determine read-count column for Escore plot; skipping plot.")
      return(NULL)
    }
  }

  row <- spreadsheet[spreadsheet$Krakenuniq.name == pathogen, ]
  threshold_escore <- if (nrow(row) > 0) row$min_escore[1] else NA
  threshold_reads <- if (nrow(row) > 0) row$min_reads[1] else NA
  
  # Get MetBrewer palette for this plot
  escore_palette <- get_met_palette("Veronese", 3)
  
  # Swap axes: Escore on Y, reads on X
  p <- ggplot(escore, aes_string(x = read_col, y = "Escore")) +
    geom_point(color = escore_palette[1], size = 6, alpha = 0.8, stroke = 2) +
    labs(title = paste("E-score vs Reads:", pathogen),
         subtitle = "Detection Metrics",
         x = "Read Count", y = "E-score") +
    beautiful_theme(18) +
    scale_x_continuous(labels = comma)
  
  if (!is.na(threshold_escore)) {
    p <- p + geom_hline(yintercept = threshold_escore, linetype = "dashed", 
                        color = escore_palette[2], linewidth = 1.5, alpha = 0.8) +
      annotate("text", x = max(escore[[read_col]], na.rm = TRUE) * 0.9, y = threshold_escore, 
               label = "E-score Threshold", hjust = -0.1, color = "black", 
               fontface = "bold", size = 5)
  }
  if (!is.na(threshold_reads)) {
    p <- p + geom_vline(xintercept = threshold_reads, linetype = "dashed", 
                        color = escore_palette[3], linewidth = 1.5, alpha = 0.8) +
      annotate("text", x = threshold_reads, y = max(escore$Escore) * 0.9, 
               label = "Read Threshold", angle = 90, vjust = -0.5, color = "black", 
               fontface = "bold", size = 5)
  }
  
  return(p)
}

# New: Guellil et al. E-value vs KrakenUniq reads plot
plot_evalue <- function() {
  if (!file.exists(escore_file)) return(NULL)
  escore <- safe_read_csv(escore_file, allow_empty = TRUE)
  if (is.null(escore) || nrow(escore) == 0) return(NULL)
  escore$taxonomy <- trimws(escore$taxonomy)
  escore <- escore[escore$taxonomy == pathogen, ]
  
  if (nrow(escore) == 0) return(NULL)
  
  # Decide which KrakenUniq read-count column to use (same logic as plot_escore)
  read_col <- NULL
  if ("reads" %in% colnames(escore)) {
    read_col <- "reads"
  } else if ("taxReads" %in% colnames(escore)) {
    read_col <- "taxReads"
  } else if ("X..reads" %in% colnames(escore)) {
    read_col <- "X..reads"
  } else if ("X.reads" %in% colnames(escore)) {
    read_col <- "X.reads"
  } else if ("X.of.reads" %in% colnames(escore)) {
    read_col <- "X.of.reads"
  } else {
    numeric_cols <- setdiff(colnames(escore), c("taxonomy"))
    numeric_cols <- numeric_cols[sapply(escore[numeric_cols], is.numeric)]
    if (length(numeric_cols) > 0) {
      read_col <- numeric_cols[1]
    } else {
      warning("Could not determine read-count column for E-value plot; skipping plot.")
      return(NULL)
    }
  }
  
  # Require Guellil_et_al_Evalue column
  if (!"Guellil_et_al_Evalue" %in% colnames(escore)) return(NULL)
  
  # Get pathogen-specific E-value threshold from spreadsheet, fallback to config default
  row <- spreadsheet[spreadsheet$Krakenuniq.name == pathogen, ]
  max_evalue <- criteria_config$guellil_evalue_threshold
  if (nrow(row) > 0) {
    for (col in c("Guellil_et_al_Evalue_threshold", "max_evalue", "evalue_threshold")) {
      if (col %in% colnames(row)) {
        val <- row[[col]][1]
        if (!is.na(val)) {
          num <- suppressWarnings(as.numeric(as.character(val)))
          if (!is.na(num)) max_evalue <- num
        }
        break
      }
    }
  }
  
  # Get read threshold from config or spreadsheet
  min_reads <- criteria_config$reads_threshold %||% 50
  if (nrow(row) > 0 && "min_reads" %in% colnames(row)) {
    val <- row$min_reads[1]
    if (!is.na(val)) {
      num <- suppressWarnings(as.numeric(as.character(val)))
      if (!is.na(num)) min_reads <- num
    }
  }
  
  # Palette
  evalue_palette <- get_met_palette("Veronese", 3)
  
  # Plot E-value (Y, log10) vs reads (X)
  p <- ggplot(escore, aes_string(x = read_col, y = "Guellil_et_al_Evalue")) +
    geom_point(color = evalue_palette[2], size = 6, alpha = 0.8, stroke = 2) +
    scale_y_log10() +
    labs(
      title = paste("Guellil et al. E-value vs Reads:", pathogen),
      subtitle = "KrakenUniq + E-value detection metric",
      x = "Read Count",
      y = "E value (log10 scale)"
    ) +
    beautiful_theme(18) +
    scale_x_continuous(labels = comma)
  
  # Horizontal line at E-value threshold (larger is better: E-value > threshold)
  if (!is.na(max_evalue) && max_evalue > 0) {
    p <- p +
      geom_hline(
        yintercept = max_evalue,
        linetype = "dashed",
        color = evalue_palette[3],
        linewidth = 1.5,
        alpha = 0.8
      )
  }
  
  # Vertical line at read threshold
  if (!is.na(min_reads) && min_reads > 0) {
    p <- p +
      geom_vline(
        xintercept = min_reads,
        linetype = "dashed",
        color = evalue_palette[3],
        linewidth = 1.5,
        alpha = 0.8
      )
  }
  
  return(p)
}

plot_ani <- function() {
  # Prefer per-read ANI distribution from BAM (created by calculate_edit_distance_r2.py)
  ani_vals <- NULL
  if (file.exists(ani_dist_file)) {
    ani_vals <- tryCatch({
      vals <- readLines(ani_dist_file)
      vals_num <- suppressWarnings(as.numeric(vals))
      vals_num[!is.na(vals_num)]
    }, error = function(e) numeric(0))
  }
  
  # Fallback: single summary ANI value from .ani.txt
  if (is.null(ani_vals) || length(ani_vals) == 0) {
  if (!file.exists(ani_file)) return(NULL)
    ani_val <- tryCatch({
      val <- readLines(ani_file, n = 1)
      suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", val)))
    }, error = function(e) NA)
  if (is.na(ani_val)) return(NULL)
    ani_vals <- c(ani_val)
  }
  
  # Get MetBrewer palette for this plot
  ani_palette <- get_met_palette("Hokusai", 4)
  
  # Compute mean ANI
  ani_mean <- mean(ani_vals, na.rm = TRUE)
  
  df <- data.frame(ANI = ani_vals)
  
  # Choose dynamic binning so dense data doesn't look like a single block
  ani_min <- min(df$ANI, na.rm = TRUE)
  ani_max <- max(df$ANI, na.rm = TRUE)
  ani_range <- max(ani_max - ani_min, 1e-6)
  
  n_reads <- nrow(df)
  # More bins when many reads; fewer when few reads
  target_bins <- if (n_reads > 5000) {
    150
  } else if (n_reads > 1000) {
    100
  } else if (n_reads > 200) {
    60
  } else {
    30
  }
  binwidth <- ani_range / target_bins
  
  p <- ggplot(df, aes(x = ANI)) +
    geom_histogram(
      binwidth = binwidth,
      fill = ani_palette[1],
      color = "black",
      linewidth = 0.6,
      alpha = 0.8,
      boundary = 0
    ) +
    geom_vline(
      xintercept = ani_mean,
      linetype = "dashed",
      color = ani_palette[3],
      linewidth = 1.5
    ) +
    labs(
      title = "ANI Distribution from BAM",
      subtitle = paste0("Mean ANI = ", round(ani_mean, 2), "%"),
      x = "ANI (%)",
      y = "Number of Reads"
    ) +
    beautiful_theme(18) +
    scale_y_continuous(labels = comma)
  
  # Sensible x-limits
  p <- p + coord_cartesian(
    xlim = c(
      max(0, floor(min(df$ANI, na.rm = TRUE) - 1)),
      min(100, ceiling(max(df$ANI, na.rm = TRUE) + 1))
    )
  )
  
  return(p)
}

# Removed entropy plot - not useful

# ------------------------- GENOME COVERAGE PLOT -------------------------
plot_genome_coverage <- function() {
  if (!file.exists(depth_file)) return(NULL)
  
  tryCatch({
    # Read depth data (one value per base position)
    depths <- scan(depth_file, what = numeric(), quiet = TRUE)
    
    if (length(depths) == 0) return(NULL)
    
    # Create data frame with positions
    genome_len <- length(depths)
    positions <- 1:genome_len
    
    # Bin data for large genomes (every 100 bases for genomes > 100kb)
    bin_size <- if (genome_len > 100000) 100 else 1
    binned_pos <- seq(1, genome_len, by = bin_size)
    binned_depth <- sapply(binned_pos, function(i) {
      end <- min(i + bin_size - 1, genome_len)
      mean(depths[i:end], na.rm = TRUE)
    })
    
    df <- data.frame(Position = binned_pos, Depth = binned_depth)
    
    # Calculate mean coverage and evenness for subtitle
    mean_cov <- mean(depths, na.rm = TRUE)
    evenness_val <- if (!is.null(summary_data) && nrow(summary_data) > 0) {
      ev <- summary_data$Evenness[1]
      if (!is.na(ev) && ev != "NA") {
        suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", ev)))
      } else NA
    } else NA
    
    subtitle_text <- paste("Mean Coverage:", round(mean_cov, 2), "X")
    if (!is.na(evenness_val)) {
      subtitle_text <- paste(subtitle_text, "| Evenness:", round(evenness_val, 2), "%")
    }
    
    # Get MetBrewer palette for this plot - Derain palette
    cov_palette <- get_met_palette("Derain", 5)
    
    p <- ggplot(df, aes(x = Position, y = Depth)) +
      geom_area(fill = cov_palette[5], alpha = 0.7, color = cov_palette[5], linewidth = 0.5) +
      geom_line(color = cov_palette[5], linewidth = 1.2) +
      labs(title = "Genome Coverage Profile",
           subtitle = subtitle_text,
           x = "Genome Position (bp)",
           y = "Coverage Depth") +
      beautiful_theme(18) +
      scale_x_continuous(labels = comma) +
      scale_y_continuous(labels = comma) +
      geom_hline(yintercept = mean_cov, linetype = "dashed", 
                color = cov_palette[1], linewidth = 1.5, alpha = 0.8) +
      annotate("text", x = max(df$Position) * 0.95, y = mean_cov, 
               label = paste("Mean:", round(mean_cov, 2)), 
               vjust = -0.5, hjust = 1, color = "black", 
               fontface = "bold", size = 5)
    
    return(p)
  }, error = function(e) {
    message("[WARNING] Could not create genome coverage plot: ", e$message)
    return(NULL)
  })
}

embed_pdf_page <- function(pdf_path, page = 1, density = 400, width = 3500) {
  if (!file.exists(pdf_path)) return(NULL)
  tryCatch({
    # Try to extract specific page and render at higher resolution
    img <- image_read_pdf(pdf_path, density = density, pages = page)
    # Resize to make larger - maintain aspect ratio
    img_info <- image_info(img)
    aspect_ratio <- img_info$height / img_info$width
    height <- as.integer(width * aspect_ratio)
    img <- image_resize(img, paste0(width, "x", height))
  rasterGrob(img, interpolate = TRUE)
  }, error = function(e) {
    message("[WARNING] Could not embed PDF page from ", pdf_path, ": ", e$message)
    return(NULL)
  })
}

# ------------------------- ADDITIONAL PLOT FUNCTIONS -------------------------
plot_coverage_metrics <- function() {
  if (is.null(summary_data) || nrow(summary_data) == 0) return(NULL)
  
  # Extract metrics with proper handling
  coverage_val <- suppressWarnings({
    cov_str <- summary_data$Coverage[1]
    if (!is.na(cov_str) && cov_str != "NA") {
      as.numeric(gsub("[^0-9\\.]", "", cov_str))
    } else NA
  })
  
  evenness_val <- suppressWarnings({
    ev <- summary_data$Evenness[1]
    if (!is.na(ev) && ev != "NA") {
      as.numeric(gsub("[^0-9\\.]", "", ev))
    } else NA
  })
  
  ani_val <- suppressWarnings({
    ani_str <- summary_data$ANI[1]
    if (!is.na(ani_str) && ani_str != "NA") {
      as.numeric(gsub("[^0-9\\.]", "", ani_str))
    } else NA
  })
  
  metrics <- data.frame(
    Metric = character(),
    Value = numeric(),
    Unit = character(),
    stringsAsFactors = FALSE
  )
  
  if (!is.na(coverage_val)) {
    metrics <- rbind(metrics, data.frame(Metric = "Coverage", Value = coverage_val, Unit = "X"))
  }
  if (!is.na(evenness_val)) {
    metrics <- rbind(metrics, data.frame(Metric = "Evenness", Value = evenness_val, Unit = "%"))
  }
  if (!is.na(ani_val)) {
    metrics <- rbind(metrics, data.frame(Metric = "ANI", Value = ani_val, Unit = "%"))
  }
  
  if (nrow(metrics) == 0) return(NULL)
  
  metrics$Label <- paste0(round(metrics$Value, 2), metrics$Unit)
  
  # Get MetBrewer palette for this plot
  cov_metrics_palette <- get_met_palette("Degas", nrow(metrics))
  
  ggplot(metrics, aes(x = reorder(Metric, Value), y = Value, fill = Metric)) +
    geom_col(alpha = 0.9, color = "black", linewidth = 1.5) +
    scale_fill_manual(values = cov_metrics_palette) +
    labs(title = "Key Quality Metrics",
         subtitle = "Coverage, Evenness, and ANI",
         x = NULL, y = "Value") +
    beautiful_theme(18) +
    theme(legend.position = "none") +
    geom_text(aes(label = Label), hjust = -0.1, size = 5, fontface = "bold") +
    coord_flip()
}

# ------------------------- DETECTION SCORE PLOT -------------------------
plot_detection_score <- function() {
  if (is.null(summary_data) || nrow(summary_data) == 0) return(NULL)
  
  # Parse Score column (format "X/Y")
  score_str <- as.character(summary_data$Score[1])
  if (is.na(score_str) || score_str == "NA" || score_str == "") return(NULL)
  
  score_parts <- strsplit(score_str, "/")[[1]]
  if (length(score_parts) != 2) return(NULL)
  
  detection_score <- suppressWarnings(as.numeric(score_parts[1]))
  max_score <- suppressWarnings(as.numeric(score_parts[2]))
  
  if (is.na(detection_score) || is.na(max_score)) return(NULL)
  
  # Parse criteria_passed (comma-separated list)
  criteria_passed_str <- as.character(summary_data$Criteria_passed[1])
  passed_criteria <- character()
  
  if (!is.na(criteria_passed_str) && 
      criteria_passed_str != "NA" && 
      criteria_passed_str != "" &&
      criteria_passed_str != "None") {
    passed_criteria <- trimws(strsplit(criteria_passed_str, ",")[[1]])
  }
  
  # Get pathogen-specific thresholds from spreadsheet
  pathogen_row <- spreadsheet[spreadsheet$Krakenuniq.name == pathogen, ]
  min_escore <- if (nrow(pathogen_row) > 0 && "min_escore" %in% colnames(pathogen_row)) {
    val <- pathogen_row$min_escore[1]
    if (!is.na(val)) as.numeric(val) else criteria_config$escore_threshold
  } else {
    criteria_config$escore_threshold
  }
  
  min_reads <- if (nrow(pathogen_row) > 0 && "min_reads" %in% colnames(pathogen_row)) {
    val <- pathogen_row$min_reads[1]
    if (!is.na(val)) as.integer(val) else criteria_config$reads_threshold
  } else {
    criteria_config$reads_threshold
  }
  
  # E-value threshold (Guellil et al.), per-pathogen or from config
  max_evalue <- criteria_config$guellil_evalue_threshold
  if (nrow(pathogen_row) > 0) {
    for (col in c("Guellil_et_al_Evalue_threshold", "max_evalue", "evalue_threshold")) {
      if (col %in% colnames(pathogen_row)) {
        val <- pathogen_row[[col]][1]
        if (!is.na(val)) {
          num <- suppressWarnings(as.numeric(as.character(val)))
          if (!is.na(num)) max_evalue <- num
        }
        break
      }
    }
  }
  
  use_evalue <- if (!is.null(criteria_config$use_evalue_for_detection)) {
    isTRUE(criteria_config$use_evalue_for_detection)
  } else {
    FALSE
  }
  
  # Check if HOPS is enabled (HOPS_score > 0 or Detected_by_HOPS is TRUE in any row)
  hops_enabled <- FALSE
  if (!is.null(summary_data) && nrow(summary_data) > 0) {
    hops_score <- suppressWarnings(as.numeric(summary_data$HOPS_score[1]))
    hops_detected <- summary_data$Detected_by_HOPS[1]
    if ((!is.na(hops_score) && hops_score > 0) || (!is.na(hops_detected) && hops_detected == TRUE)) {
      hops_enabled <- TRUE
    }
  }
  
  # Get all possible criteria from config and summary data (with pathogen-specific thresholds)
  # First criterion: either E-score or E-value depending on config
  if (use_evalue) {
    all_criteria <- c(
      paste0("E-value > ", max_evalue),
      paste0("Reads >= ", min_reads)
    )
  } else {
    all_criteria <- c(
      paste0("E-Score >= ", min_escore),
      paste0("Reads >= ", min_reads)
    )
  }
  
  # Only add HOPS criteria if HOPS is enabled
  if (hops_enabled) {
    all_criteria <- c(all_criteria,
      "HOPS Detection",
      "HOPS Edit Distance",
      "HOPS Damage"
    )
  }
  
  # Get read mapping ratio threshold from config
  read_ratio_thresh <- if (!is.null(criteria_config$read_mapping_ratio_threshold)) {
    criteria_config$read_mapping_ratio_threshold
  } else {
    0.5  # Default threshold
  }
  
  # Add remaining criteria
  edit_damage_available <- file.exists(edit_r2_damaged_file) && file.exists(edit_r2_no_damage_file)
  edit_decay_thresh <- if (!is.null(criteria_config$edit_distance_decay_quality_threshold)) {
    criteria_config$edit_distance_decay_quality_threshold
  } else {
    criteria_config$edit_distance_logr2_threshold %||% 0.65
  }
  edit_no_damage_tile <- if (!is.null(criteria_config$edit_distance_no_damage_decay_quality_min)) {
    paste0("Edit Distance No-Damage Decay Quality >= ", criteria_config$edit_distance_no_damage_decay_quality_min)
  } else if (!is.null(criteria_config$edit_distance_no_damage_decay_quality_max)) {
    paste0("Edit Distance No-Damage Decay Quality <= ", criteria_config$edit_distance_no_damage_decay_quality_max)
  } else {
    paste0("Edit Distance No-Damage Decay Quality >= ", 0.55)
  }

  edit_criteria_tiles <- if (edit_damage_available) {
    c(
      paste0("Edit Distance Damage Decay Quality >= ", edit_decay_thresh),
      edit_no_damage_tile
    )
  } else {
    c(
      paste0("Edit Distance Decay Quality >= ", edit_decay_thresh)
    )
  }

  all_criteria <- c(
    all_criteria,
    paste0("ANI > ", criteria_config$ani_threshold, "%"),
    paste0("5' C>T Damage >= ", criteria_config$damage_5p_ct_threshold),
    paste0("Breadth Ratio >= ", criteria_config$breadth_ratio_threshold),
    # Entropy threshold label will be dynamically determined based on virus/bacteria
    # We'll match it from the passed criteria in the summary CSV
    "Entropy",
    edit_criteria_tiles,
    "Genus Ranking = 1",
    paste0("Read Mapping Ratio >= ", read_ratio_thresh)
  )
  
  # Create criteria list with pass/fail status
  criteria_list <- character()
  passed_list <- logical()
  
  for (criterion in all_criteria) {
    # Check if this criterion is in the passed list (exact match or contains the criterion name)
    # Handle both exact matches and partial matches (e.g., "E-Score >= 5" matches "E-Score >= 5")
    # Special handling for "Entropy" - match any entropy criterion regardless of threshold or (virus)/(bacteria) suffix
    if (criterion == "Entropy") {
      # Match any entropy criterion (e.g., "Entropy (1000bp) >= 0.9 (bacteria)" or "Entropy (1000bp) >= 0.7 (virus)")
      is_passed <- any(grepl("Entropy.*>=", passed_criteria, ignore.case = TRUE))
    } else if (startsWith(criterion, "Edit Distance Damage Decay Quality")) {
      # Do NOT match the numeric threshold exactly: config thresholds and summary labels
      # can differ in formatting (0.65 vs 0.6500) or source (config.yaml vs spreadsheet).
      # Treat as passed if any passed criterion mentions this edit-distance metric.
      is_passed <- any(grepl("^Edit Distance Damage Decay Quality\\s*>=", passed_criteria)) ||
        any(grepl("^Edit Distance Decay Quality\\s*>=", passed_criteria))
    } else if (startsWith(criterion, "Edit Distance No-Damage Decay Quality")) {
      # Same idea: match by prefix rather than exact numeric.
      is_passed <- any(grepl("^Edit Distance No-Damage Decay Quality\\s*>=", passed_criteria)) ||
        any(grepl("^Edit Distance No-Damage Decay Quality\\s*<=", passed_criteria)) ||
        any(grepl("^Edit Distance Decay Quality\\s*>=", passed_criteria))
    } else {
      is_passed <- any(grepl(criterion, passed_criteria, fixed = TRUE))
    }
    criteria_list <- c(criteria_list, criterion)
    passed_list <- c(passed_list, is_passed)
  }
  
  if (length(criteria_list) == 0) {
    return(NULL)  # No applicable criteria
  }
  
  score_breakdown <- data.frame(
    Criterion = criteria_list,
    Passed = passed_list,
    Value = ifelse(passed_list, 1, 0),
    stringsAsFactors = FALSE
  )
  
  # Binary heatmap for criteria with more distinct colors
  score_breakdown$Criterion <- factor(score_breakdown$Criterion, 
                                       levels = rev(score_breakdown$Criterion))
  
  # Get MetBrewer palette for detection score heatmap - Hokusai1
  score_palette <- get_met_palette("Hokusai1", 6)
  
  # Use Hokusai1 5th color (green) for pass, Hokusai1 2nd color for fail
  p2 <- ggplot(score_breakdown, aes(x = 1, y = Criterion, fill = Value)) +
    geom_tile(color = "black", linewidth = 1.5, width = 0.5) +
    scale_fill_gradient2(low = score_palette[2], mid = "#FFC107", high = score_palette[5], 
                        midpoint = 0.5, limits = c(0, 1), guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = "Pathogen Detection Score",
         subtitle = paste("Score:", detection_score, "/", max_score, 
                         paste0("(", round(detection_score/max_score*100, 1), "%)")),
         x = NULL, y = NULL) +
    beautiful_theme(18) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(10, 20, 10, 20),
          plot.title = element_text(size = 20, face = "bold"),
          plot.subtitle = element_text(size = 16)) +
    geom_text(aes(label = ifelse(Value == 1, "PASS", "FAIL")), 
              size = 6, fontface = "bold", 
              color = "black")
  
  # Return only the heatmap (removed the bar chart)
  return(list(p2))
}

# ------------------------- RELATIVE ENTROPY AND BREADTH RATIO PLOT (Nature paper style) -------------------------
plot_entropy_breadth_ratio <- function() {
  if (is.null(summary_data) || nrow(summary_data) == 0) return(NULL)
  
  # Get relative entropy values (pathopipe method: 100bp and 1000bp windows)
  entropy_100_val <- suppressWarnings({
    ent <- summary_data$Relative_entropy_100bp[1]
    if (is.numeric(ent) && !is.na(ent)) ent else NA
  })
  entropy_1000_val <- suppressWarnings({
    ent <- summary_data$Relative_entropy_1000bp[1]
    if (is.numeric(ent) && !is.na(ent)) ent else NA
  })
  
  # Fallback to legacy file if new columns not available
  if (is.na(entropy_100_val) || is.na(entropy_1000_val)) {
    entropy_legacy_file_100 <- gsub("\\.mean_entropy\\.txt$", ".entropy_100bp.txt", entropy_mean_file)
    entropy_legacy_file_1000 <- gsub("\\.mean_entropy\\.txt$", ".entropy_1000bp.txt", entropy_mean_file)
    
    if (is.na(entropy_100_val) && file.exists(entropy_legacy_file_100)) {
      entropy_100_val <- suppressWarnings(as.numeric(read_metric_file(entropy_legacy_file_100)))
    }
    if (is.na(entropy_1000_val) && file.exists(entropy_legacy_file_1000)) {
      entropy_1000_val <- suppressWarnings(as.numeric(read_metric_file(entropy_legacy_file_1000)))
    }
    # Final fallback to legacy mean_entropy.txt (assumes 1000bp value)
    if (is.na(entropy_1000_val)) {
      entropy_1000_val <- suppressWarnings(as.numeric(read_metric_file(entropy_mean_file)))
    }
  }
  
  # Get breadth ratio
  breadth_ratio_val <- suppressWarnings({
    br <- summary_data$Breadth_ratio[1]
    if (is.numeric(br) && !is.na(br)) br else NA
  })
  if (is.na(breadth_ratio_val)) {
    breadth_ratio_val <- suppressWarnings(read_metric_file(breadth_ratio_file))
  }
  
  # Determine if this is a virus for threshold selection
  # According to Nature paper: ALL viruses use 0.7 threshold (not just short ones)
  # Detection: check if pathogen name contains "virus" (case-insensitive)
  is_virus <- FALSE
  if (!is.null(summary_data) && nrow(summary_data) > 0) {
    pathogen_name_check <- as.character(summary_data$Pathogen[1])
    if (grepl("virus", pathogen_name_check, ignore.case = TRUE)) {
      is_virus <- TRUE
    }
  }
  # Also check BAM file to get reference length as additional check
  # (short genomes < 10kb are also treated as viruses)
  if (file.exists(bam_file)) {
    tryCatch({
      if (requireNamespace("Rsamtools", quietly = TRUE)) {
        bam_handle <- Rsamtools::BamFile(bam_file)
        ref_info <- Rsamtools::seqinfo(bam_handle)
        if (length(ref_info@seqlengths) > 0) {
          ref_len <- ref_info@seqlengths[1]
          if (!is.na(ref_len) && ref_len < 10000) {
            is_virus <- TRUE
          }
        }
      }
    }, error = function(e) {
      # If we can't determine, rely on name-based detection
    })
  }
  
  # Create data frame with all metrics
  metric_names <- c()
  metric_values <- c()
  metric_available <- c()
  
  if (!is.na(entropy_100_val)) {
    metric_names <- c(metric_names, "Entropy (100bp)")
    metric_values <- c(metric_values, entropy_100_val)
    metric_available <- c(metric_available, TRUE)
  }
  if (!is.na(entropy_1000_val)) {
    metric_names <- c(metric_names, "Entropy (1000bp)")
    metric_values <- c(metric_values, entropy_1000_val)
    metric_available <- c(metric_available, TRUE)
  }
  if (!is.na(breadth_ratio_val)) {
    metric_names <- c(metric_names, "Breadth Ratio")
    metric_values <- c(metric_values, breadth_ratio_val)
    metric_available <- c(metric_available, TRUE)
  }
  
  if (length(metric_names) == 0) return(NULL)
  
  df <- data.frame(
    Metric = metric_names,
    Value = metric_values,
    Available = metric_available
  )
  df <- df[df$Available, ]
  if (nrow(df) == 0) return(NULL)
  
  # Get thresholds from config
  # Entropy threshold: 0.9 for bacteria, 0.7 for viruses (all viruses, not just short ones)
  if (is_virus) {
    entropy_thresh <- if (!is.null(criteria_config$entropy_threshold_virus)) {
      criteria_config$entropy_threshold_virus
    } else {
      0.7
    }
  } else {
    entropy_thresh <- if (!is.null(criteria_config$entropy_threshold)) {
      criteria_config$entropy_threshold
    } else {
      0.9
    }
  }
  breadth_thresh <- criteria_config$breadth_ratio_threshold %||% 0.8
  
  # Get MetBrewer palette for this plot - Morgenstern palette
  eb_palette <- get_met_palette("Morgenstern", 3)
  
  # Map colors: entropy 100bp, entropy 1000bp, breadth ratio
  fill_colors <- c(
    "Entropy (100bp)" = eb_palette[1],
    "Entropy (1000bp)" = eb_palette[2],
    "Breadth Ratio" = eb_palette[3]
  )
  
  # Calculate mean entropy for display
  mean_entropy <- NA
  if (!is.na(entropy_100_val) && !is.na(entropy_1000_val)) {
    mean_entropy <- mean(c(entropy_100_val, entropy_1000_val))
  } else if (!is.na(entropy_1000_val)) {
    mean_entropy <- entropy_1000_val
  } else if (!is.na(entropy_100_val)) {
    mean_entropy <- entropy_100_val
  }
  
  # Create separate data frames for entropy and breadth ratio for faceting
  entropy_df <- data.frame(
    Metric = character(),
    Value = numeric(),
    Type = character(),
    stringsAsFactors = FALSE
  )
  
  if (!is.na(entropy_100_val)) {
    entropy_df <- rbind(entropy_df, data.frame(
      Metric = "Entropy (100bp)",
      Value = entropy_100_val,
      Type = "Entropy"
    ))
  }
  if (!is.na(entropy_1000_val)) {
    entropy_df <- rbind(entropy_df, data.frame(
      Metric = "Entropy (1000bp)",
      Value = entropy_1000_val,
      Type = "Entropy"
    ))
  }
  
  breadth_df <- NULL
  if (!is.na(breadth_ratio_val)) {
    breadth_df <- data.frame(
      Metric = "Breadth Ratio",
      Value = breadth_ratio_val,
      Type = "Breadth Ratio"
    )
  }
  
  # Combine for faceting
  if (nrow(entropy_df) > 0 && !is.null(breadth_df)) {
    combined_df <- rbind(entropy_df, breadth_df)
  } else if (nrow(entropy_df) > 0) {
    combined_df <- entropy_df
  } else if (!is.null(breadth_df)) {
    combined_df <- breadth_df
  } else {
    return(NULL)
  }
  
  # Create facet labels
  combined_df$Facet <- ifelse(combined_df$Type == "Entropy", "Relative Entropy", "Breadth Ratio")
  combined_df$Metric <- factor(combined_df$Metric, levels = unique(combined_df$Metric))
  
  # Create faceted plot
  p <- ggplot(combined_df, aes(x = Metric, y = Value, fill = Metric)) +
    geom_col(alpha = 0.9, color = "black", linewidth = 2, width = 0.6) +
    scale_fill_manual(values = fill_colors) +
    facet_wrap(~ Facet, scales = "free_x", ncol = 2) +
    labs(title = "Quality Metrics: Relative Entropy & Breadth Ratio",
         subtitle = if (!is.na(mean_entropy)) {
           paste0("Mean Entropy: ", round(mean_entropy, 4))
         } else {
           ""
         },
         x = NULL, y = "Value") +
    beautiful_theme(18) +
    theme(legend.position = "none") +
    geom_text(aes(label = round(Value, 4)), vjust = -0.5, size = 6, fontface = "bold")
  
  # Create threshold data frames for each facet
  threshold_df <- data.frame(
    Facet = character(),
    Threshold = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (nrow(entropy_df) > 0) {
    threshold_df <- rbind(threshold_df, data.frame(
      Facet = "Relative Entropy",
      Threshold = entropy_thresh
    ))
  }
  if (!is.null(breadth_df)) {
    threshold_df <- rbind(threshold_df, data.frame(
      Facet = "Breadth Ratio",
      Threshold = breadth_thresh
    ))
  }
  
  # Add threshold lines using geom_hline with facet data
  if (nrow(threshold_df) > 0) {
    p <- p +
      geom_hline(
        data = threshold_df,
        aes(yintercept = Threshold),
        linetype = "dashed",
        color = "black",
        linewidth = 1.5,
        alpha = 0.8
      )
  }
  
  # Adjust y-axis - calculate max for both facets
  max_val_entropy <- if (nrow(entropy_df) > 0) {
    max(c(entropy_df$Value, entropy_thresh, 1.0), na.rm = TRUE)
  } else 1.0
  max_val_breadth <- if (!is.null(breadth_df)) {
    max(c(breadth_df$Value, breadth_thresh, 1.0), na.rm = TRUE)
  } else 1.0
  
  max_val <- max(max_val_entropy, max_val_breadth, na.rm = TRUE) * 1.2
  p <- p + ylim(0, max_val)
  
  return(p)
}

plot_damage_metrics <- function() {
  if (is.null(summary_data) || nrow(summary_data) == 0) return(NULL)
  
  damage_5p <- suppressWarnings(as.numeric(summary_data$Damage_5p_CtoT[1]))
  entropy <- suppressWarnings(as.numeric(read_metric_file(entropy_mean_file)))
  
  if (is.na(damage_5p) && is.na(entropy)) return(NULL)
  
  metrics <- data.frame(
    Metric = c("5' C>T Damage", "Mean Entropy"),
    Value = c(ifelse(is.na(damage_5p), 0, damage_5p), ifelse(is.na(entropy), 0, entropy))
  )
  metrics <- metrics[metrics$Value > 0, ]
  if (nrow(metrics) == 0) return(NULL)
  
  # Get MetBrewer palette for this plot
  damage_palette <- get_met_palette("Monet", 2)
  
  ggplot(metrics, aes(x = Metric, y = Value, fill = Metric)) +
    geom_col(alpha = 0.9, color = "black", linewidth = 1.5, width = 0.6) +
    scale_fill_manual(values = damage_palette) +
    labs(title = "Damage and Entropy Metrics",
         subtitle = "Ancient DNA Characteristics",
         x = NULL, y = "Value") +
    beautiful_theme(18) +
    theme(legend.position = "none") +
    geom_text(aes(label = round(Value, 4)), vjust = -0.5, size = 6, fontface = "bold")
}

# ------------------------- EDIT DISTANCE PLOT (Occurrences vs Mismatches 0-5) -------------------------
plot_edit_distance_distribution <- function() {
  make_single_plot <- function(dist_file, r2_file, subset_label) {
    if (!file.exists(dist_file)) return(NULL)

    ed_dist <- tryCatch({
      read.table(dist_file, sep = "\t", header = FALSE, col.names = c("Mismatches", "Occurrences"))
    }, error = function(e) NULL)

    if (is.null(ed_dist) || nrow(ed_dist) == 0) return(NULL)

    # Filter to first 6 (0-5)
    ed_dist <- ed_dist[ed_dist$Mismatches <= 5, ]

    # Read Decay Quality Score and regression metrics from r2 file
    decay_quality_score <- NA
    r2_linear <- NA
    r2_log <- NA
    edit_r2_val <- NA

    if (!is.null(r2_file) && file.exists(r2_file)) {
      lines <- readLines(r2_file)
      if (length(lines) >= 1) decay_quality_score <- suppressWarnings(as.numeric(lines[1]))

      for (line in lines) {
        line <- trimws(line)
        if (grepl("^Legacy_metric=", line)) {
          edit_r2_val <- suppressWarnings(as.numeric(gsub("Legacy_metric=", "", line)))
        } else if (grepl("^R2_linear=", line)) {
          r2_linear <- suppressWarnings(as.numeric(gsub("R2_linear=", "", line)))
        } else if (grepl("^R2=", line)) {
          r2_linear <- suppressWarnings(as.numeric(gsub("R2=", "", line)))
        } else if (grepl("^R2_log=", line)) {
          r2_log <- suppressWarnings(as.numeric(gsub("R2_log=", "", line)))
        } else if (grepl("^Decay_quality_score=", line)) {
          if (is.na(decay_quality_score)) {
            decay_quality_score <- suppressWarnings(as.numeric(gsub("Decay_quality_score=", "", line)))
          }
        }
      }
    }

    # Two visual styles: damage (red) vs no-damage (blue)
    is_damage <- tolower(subset_label) == "damage"
    bar_fill <- if (is_damage) "#E15261" else "#4477AA"
    lin_color <- if (is_damage) "#b83a3f" else "#2a63a6"
    log_color <- if (is_damage) "#f18a8f" else "#6baed6"

    # Fit regression lines from distribution
    regression_data <- NULL
    ed_dist_nonzero <- ed_dist[ed_dist$Occurrences > 0, ]
    if (nrow(ed_dist_nonzero) >= 2) {
      lm_lin <- lm(Occurrences ~ Mismatches, data = ed_dist_nonzero)
      lm_log <- lm(log1p(Occurrences) ~ Mismatches, data = ed_dist_nonzero)

      regression_data <- data.frame(Mismatches = seq(0, 5, by = 0.1))
      regression_data$Occ_linear <- predict(lm_lin, newdata = regression_data)
      regression_data$Occ_log <- exp(predict(lm_log, newdata = regression_data)) - 1

      regression_data$Occ_linear <- pmax(0, regression_data$Occ_linear)
      regression_data$Occ_log <- pmax(0, regression_data$Occ_log)
      max_occ_val <- max(ed_dist$Occurrences)
      regression_data$Occ_linear <- pmin(regression_data$Occ_linear, max_occ_val * 2)
      regression_data$Occ_log <- pmin(regression_data$Occ_log, max_occ_val * 2)
    }

    max_occ <- max(ed_dist$Occurrences)
    if (!is.null(regression_data) && nrow(regression_data) > 0) {
      max_regression <- max(regression_data$Occ_linear, regression_data$Occ_log, na.rm = TRUE)
      max_occ <- max(max_occ, max_regression, na.rm = TRUE)
    }
    if (max_occ == 0) max_occ <- 1
    y_max <- max_occ * 1.2

    p <- ggplot(ed_dist, aes(x = Mismatches, y = Occurrences)) +
      geom_col(fill = bar_fill, alpha = 0.9, color = "black", linewidth = 1.5, width = 0.7) +
      labs(
        title = if (is_damage) "Edit distance (Damage)" else "Edit distance",
        subtitle = if (!is.na(decay_quality_score)) {
          paste0("Decay Quality Score = ", round(decay_quality_score, 4))
        } else {
          ""
        },
        x = "Number of Mismatches", y = "Number of Occurrences"
      ) +
      beautiful_theme(20) +
      scale_x_continuous(breaks = 0:5, labels = 0:5) +
      scale_y_continuous(labels = comma, limits = c(0, y_max)) +
      theme(legend.position = "none")

    # Add regression lines
    if (!is.null(regression_data) && nrow(regression_data) > 0) {
      if (any(!is.na(regression_data$Occ_linear))) {
        lin_clean <- regression_data[!is.na(regression_data$Occ_linear), ]
        if (nrow(lin_clean) > 0) {
          p <- p + geom_line(
            data = lin_clean,
            aes(x = Mismatches, y = Occ_linear),
            linetype = "dashed",
            color = lin_color,
            linewidth = 1.3,
            alpha = 0.8
          )
        }
      }
      if (any(!is.na(regression_data$Occ_log))) {
        log_clean <- regression_data[!is.na(regression_data$Occ_log), ]
        if (nrow(log_clean) > 0) {
          p <- p + geom_line(
            data = log_clean,
            aes(x = Mismatches, y = Occ_log),
            linetype = "dotted",
            color = log_color,
            linewidth = 1.3,
            alpha = 0.8
          )
        }
      }
    }

    return(p)
  }

  # Prefer damage/no-damage metrics when present (HOPS-like).
  if (file.exists(edit_r2_damaged_dist_file) && file.exists(edit_r2_damaged_file)) {
    p_damage <- make_single_plot(edit_r2_damaged_dist_file, edit_r2_damaged_file, "damage")
    p_no_damage <- NULL
    if (file.exists(edit_r2_no_damage_dist_file) && file.exists(edit_r2_no_damage_file)) {
      p_no_damage <- make_single_plot(edit_r2_no_damage_dist_file, edit_r2_no_damage_file, "no-damage")
    }
    # Prefer: left = no-damage, right = damage
    if (!is.null(p_no_damage)) {
      return(arrangeGrob(
        grobs = list(p_no_damage, p_damage),
        ncol = 2,
        padding = unit(0.5, "cm")
      ))
    }
    return(p_damage)
  }

  # Fallback to legacy single edit-distance plot
  return(make_single_plot(edit_r2_dist_file, edit_r2_file, "all reads"))
}

# ------------------------- GENUS RANKING BAR PLOT -------------------------
plot_genus_ranking <- function() {
  if (is.null(summary_data) || nrow(summary_data) == 0) return(NULL)
  
  # Get genus ranking from summary or file (handle "#5" format)
  genus_rank_str <- summary_data$Genus_ranking[1]
  genus_rank <- NA
  
  if (!is.na(genus_rank_str) && genus_rank_str != "NA" && genus_rank_str != "error") {
    # Handle "#5" format
    if (is.character(genus_rank_str) && grepl("^#", genus_rank_str)) {
      genus_rank <- suppressWarnings(as.numeric(gsub("#", "", genus_rank_str)))
    } else {
      genus_rank <- suppressWarnings(as.numeric(genus_rank_str))
    }
  }
  
  if (is.na(genus_rank) && file.exists(genus_ranking_file)) {
    genus_rank <- suppressWarnings(as.numeric(read_metric_file(genus_ranking_file)))
  }
  
  if (is.na(genus_rank)) return(NULL)
  
  # Read KrakenUniq report to get all species in the genus with their read counts
  if (!file.exists(kraken_report_file)) return(NULL)
  
  # Parse KrakenUniq report
  kraken_lines <- readLines(kraken_report_file)
  header_idx <- which(grepl("^%", kraken_lines))[1]
  if (is.na(header_idx)) return(NULL)
  
  # Read the report data
  header <- strsplit(kraken_lines[header_idx], "\t")[[1]]
  header <- gsub("^%", "", header)
  header <- trimws(header)
  
  # Parse data lines
  data_lines <- kraken_lines[(header_idx + 1):length(kraken_lines)]
  if (length(data_lines) == 0) return(NULL)
  
  # Create a simple data frame from the report
  report_data <- tryCatch({
    read.table(text = paste(data_lines, collapse = "\n"), 
                sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                fill = TRUE, quote = "")
  }, error = function(e) NULL)
  
  if (is.null(report_data) || nrow(report_data) == 0) return(NULL)
  
  # Set column names (adjust based on actual header)
  n_cols <- min(length(header), ncol(report_data))
  colnames(report_data)[1:n_cols] <- header[1:n_cols]
  
  # Find taxName or taxonomy column
  tax_col <- NULL
  for (col in c("taxName", "taxonomy", "taxName", "name")) {
    if (col %in% colnames(report_data)) {
      tax_col <- col
      break
    }
  }
  if (is.null(tax_col)) return(NULL)
  
  # Find reads column
  reads_col <- NULL
  for (col in c("reads", "taxReads", "Reads")) {
    if (col %in% colnames(report_data)) {
      reads_col <- col
      break
    }
  }
  if (is.null(reads_col)) return(NULL)
  
  # Find rank column
  rank_col <- NULL
  for (col in c("rank", "Rank")) {
    if (col %in% colnames(report_data)) {
      rank_col <- col
      break
    }
  }
  if (is.null(rank_col)) return(NULL)
  
  # Get pathogen name from summary
  pathogen_name <- as.character(summary_data$Pathogen[1])
  
  # Try to read genus species from the ranking file (new format includes all species)
  genus_species <- NULL
  if (file.exists(genus_ranking_file)) {
    ranking_lines <- readLines(genus_ranking_file)
    if (length(ranking_lines) > 1) {
      # New format: first line is ranking, rest are species (taxID|taxonomy|reads)
      genus_species_list <- list()
      for (i in 2:length(ranking_lines)) {
        parts <- strsplit(ranking_lines[i], "\\|")[[1]]
        if (length(parts) >= 3) {
          genus_species_list[[length(genus_species_list) + 1]] <- data.frame(
            taxID = as.numeric(parts[1]),
            taxonomy = parts[2],
            reads = as.numeric(parts[3]),
            stringsAsFactors = FALSE
          )
        }
      }
      if (length(genus_species_list) > 0) {
        genus_species <- do.call(rbind, genus_species_list)
      }
    }
  }
  
  # Fallback: parse from KrakenUniq report using hierarchy (like Python script)
  if (is.null(genus_species) || nrow(genus_species) == 0) {
    # Fallback: create simple display with just the ranking
    ranking_text <- paste0("Genus Ranking: #", genus_rank)
    if (genus_rank == 1) {
      ranking_text <- paste0(ranking_text, " (Rank 1 - PASS)")
    } else {
      ranking_text <- paste0(ranking_text, " (FAIL)")
    }
    
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = ranking_text, 
               size = 14, fontface = "bold", hjust = 0.5, vjust = 0.5,
               color = "black") +
      xlim(0, 1) + ylim(0, 1) +
      theme_void() +
      labs(title = "Genus Ranking Position",
           subtitle = "Position of pathogen within its genus (by read count). Rank #1 = PASS")
    return(p)
  }
  
  # Convert reads to numeric (if not already)
  if (!"reads_num" %in% colnames(genus_species)) {
    if ("reads" %in% colnames(genus_species)) {
      genus_species$reads_num <- genus_species$reads
    } else if (reads_col %in% colnames(genus_species)) {
      genus_species$reads_num <- genus_species[[reads_col]]
    } else {
      return(NULL)
    }
  }
  genus_species$reads_num <- suppressWarnings(as.numeric(genus_species$reads_num))
  genus_species <- genus_species[!is.na(genus_species$reads_num), ]
  
  if (nrow(genus_species) == 0) return(NULL)
  
  # Sort by reads (descending - highest first)
  genus_species <- genus_species[order(-genus_species$reads_num), ]
  
  # Create species labels (truncate long names)
  if (!"taxonomy" %in% colnames(genus_species)) {
    if (tax_col %in% colnames(genus_species)) {
      genus_species$taxonomy <- genus_species[[tax_col]]
    } else {
      return(NULL)
    }
  }
  genus_species$species_label <- substr(genus_species$taxonomy, 1, 40)
  genus_species$species_label <- trimws(genus_species$species_label)
  
  # FIXED: Handle duplicate labels by making them unique
  if (any(duplicated(genus_species$species_label))) {
    genus_species$species_label <- make.unique(genus_species$species_label, sep = "_")
  }
  
  # Identify target pathogen (highlight it)
  genus_species$is_target <- grepl(pathogen_name, genus_species$taxonomy, ignore.case = TRUE)
  
  # Get MetBrewer palette - Cassat2 palette
  rank_palette <- get_met_palette("Cassatt2", 4)
  # Use different colors: target pathogen gets palette[1] (highlight color)
  # Other species get palette[4] (default color)
  
  # Create fill color column based on is_target
  genus_species$fill_color <- ifelse(genus_species$is_target, rank_palette[1], rank_palette[4])
  
  # Create factor for ordering - highest at top after coord_flip()
  unique_labels <- unique(genus_species$species_label)
  genus_species$species_label_ordered <- factor(genus_species$species_label, 
                                                levels = rev(unique_labels))
  
  # Create custom labels with markdown formatting for bold target pathogen
  tryCatch({
    if (requireNamespace("ggtext", quietly = TRUE)) {
      library(ggtext)
      genus_species$label_formatted <- ifelse(
        genus_species$is_target,
        paste0("**", genus_species$species_label, "**"),
        genus_species$species_label
      )
      
      label_map <- setNames(genus_species$label_formatted, as.character(genus_species$species_label_ordered))
      
      p <- ggplot(genus_species, aes(x = species_label_ordered, y = reads_num)) +
        geom_col(aes(fill = fill_color), alpha = 0.9, color = "black", linewidth = 1.5) +
        scale_fill_identity(guide = "none") +
        scale_x_discrete(labels = label_map) +
        labs(title = paste0("Genus Ranking: #", genus_rank, if (genus_rank == 1) " (PASS)" else " (FAIL)"),
             x = "Species", y = "Krakenuniq reads") +
        beautiful_theme(20) +
        theme(axis.text.y = element_markdown(angle = 0, hjust = 1, size = 12)) +
        scale_y_continuous(labels = comma) +
        coord_flip()
    } else {
      stop("ggtext not available")
    }
  }, error = function(e) {
    # Fallback: use plain text (target pathogen highlighted by different fill color)
    p <<- ggplot(genus_species, aes(x = species_label_ordered, y = reads_num)) +
      geom_col(aes(fill = fill_color), alpha = 0.9, color = "black", linewidth = 1.5) +
      scale_fill_identity(guide = "none") +
      labs(title = paste0("Genus Ranking: #", genus_rank, if (genus_rank == 1) " (PASS)" else " (FAIL)"),
           x = "Species", y = "Krakenuniq reads") +
      beautiful_theme(20) +
      theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12)) +
      scale_y_continuous(labels = comma) +
      coord_flip()
  })
  
  return(p)
}

# ------------------------- CREATE PLOTS -------------------------
title_page <- create_title_page()
plot_escore_obj <- plot_escore()
plot_evalue_obj <- plot_evalue()
plot_ani_obj <- plot_ani()
plot_genome_cov_obj <- plot_genome_coverage()
plot_coverage_metrics_obj <- plot_coverage_metrics()
plot_damage_metrics_obj <- plot_damage_metrics()
plot_detection_score_obj <- plot_detection_score()
# AdnaPlotter removed - too time and memory consuming
damage_mis_grob <- embed_pdf_page(damage_pdfs$mis, density = 400, width = 3500)
damage_edit_grob <- embed_pdf_page(damage_pdfs$edit, density = 400, width = 3500)
damage_length_grob <- embed_pdf_page(damage_pdfs$length, density = 400, width = 3500)
hops_grob <- if (!is.null(summary_pdf) && file.exists(summary_pdf)) embed_pdf_page(summary_pdf, density = 400, width = 3500) else NULL
plot_entropy_breadth_obj <- plot_entropy_breadth_ratio()
plot_edit_dist_obj <- plot_edit_distance_distribution()
plot_genus_rank_obj <- plot_genus_ranking()

# ------------------------- CREATE MULTI-PAGE PDF -------------------------
all_pages <- list()

# Page 1: Title page
if (!is.null(title_page)) {
  all_pages[[length(all_pages) + 1]] <- title_page
}

# Page 2: Detection Score (heatmap only with score in subtitle)
if (!is.null(plot_detection_score_obj) && is.list(plot_detection_score_obj)) {
  page2_layout <- arrangeGrob(grobs = plot_detection_score_obj, ncol = 1,
                              padding = unit(1, "cm"))
  all_pages[[length(all_pages) + 1]] <- page2_layout
}

# Page 3: Detection metrics (E-score or E-value based on config, HOPS if available, ANI)
page3_grobs <- list()
# Only plot E-value when E-value mode is activated, otherwise use Escore
if (isTRUE(criteria_config$use_evalue_for_detection)) {
  if (!is.null(plot_evalue_obj)) page3_grobs[[length(page3_grobs) + 1]] <- plot_evalue_obj
} else {
  if (!is.null(plot_escore_obj)) page3_grobs[[length(page3_grobs) + 1]] <- plot_escore_obj
}
if (!is.null(hops_grob)) page3_grobs[[length(page3_grobs) + 1]] <- hops_grob
if (!is.null(plot_ani_obj)) page3_grobs[[length(page3_grobs) + 1]] <- plot_ani_obj

if (length(page3_grobs) > 0) {
  page3_layout <- arrangeGrob(grobs = page3_grobs, ncol = min(2, length(page3_grobs)), 
                              padding = unit(1, "cm"),
                              top = textGrob("Detection Metrics", 
                                           gp = gpar(fontsize = 20, fontface = "bold", col = "#2c3e50"),
                                           just = "top"))
  all_pages[[length(all_pages) + 1]] <- page3_layout
}

# Page 4: Genus Ranking (full page) – directly after E-score + ANI
if (!is.null(plot_genus_rank_obj)) {
  page4_layout <- arrangeGrob(
    grobs = list(plot_genus_rank_obj), 
    padding = unit(1, "cm"),
    top = textGrob(
      "Genus Ranking", 
      gp = gpar(fontsize = 20, fontface = "bold", col = "#2c3e50"),
      just = "top"
    )
  )
  all_pages[[length(all_pages) + 1]] <- page4_layout
}

# Page 5: Genome coverage plot (full page for better visibility)
if (!is.null(plot_genome_cov_obj)) {
  page5_layout <- arrangeGrob(
    grobs = list(plot_genome_cov_obj), 
    padding = unit(1, "cm"),
    top = textGrob(
      "Genome Coverage Analysis", 
      gp = gpar(fontsize = 20, fontface = "bold", col = "#2c3e50"),
      just = "top"
    )
  )
  all_pages[[length(all_pages) + 1]] <- page5_layout
}

# Page 6: Quality metrics (Relative Entropy & Breadth Ratio as per Nature paper)
if (!is.null(plot_entropy_breadth_obj)) {
  page6_layout <- arrangeGrob(
    grobs = list(plot_entropy_breadth_obj), 
    padding = unit(1, "cm"),
    top = textGrob(
      "Quality Metrics: Relative Entropy & Breadth Ratio", 
      gp = gpar(fontsize = 20, fontface = "bold", col = "#2c3e50"),
      just = "top"
    )
  )
  all_pages[[length(all_pages) + 1]] <- page6_layout
}

# Page 7: Edit Distance Distribution (full page)
if (!is.null(plot_edit_dist_obj)) {
  page7_layout <- arrangeGrob(
    grobs = list(plot_edit_dist_obj), 
    padding = unit(1, "cm"),
    top = textGrob(
      "Edit Distance Distribution", 
      gp = gpar(fontsize = 20, fontface = "bold", col = "#2c3e50"),
      just = "top"
    )
  )
  all_pages[[length(all_pages) + 1]] <- page7_layout
}

# Page 8+: Damage profiler plots - split across 2 pages for better visibility
damage_grobs <- list()
if (!is.null(damage_mis_grob)) damage_grobs[[length(damage_grobs) + 1]] <- damage_mis_grob
if (!is.null(damage_edit_grob)) damage_grobs[[length(damage_grobs) + 1]] <- damage_edit_grob
if (!is.null(damage_length_grob)) damage_grobs[[length(damage_grobs) + 1]] <- damage_length_grob

if (length(damage_grobs) > 0) {
  # Split into 2 pages: first 2 plots on one page, remaining on the next
  if (length(damage_grobs) >= 2) {
    # First damage page
    damage_page1 <- arrangeGrob(
      grobs = damage_grobs[1:2], 
      ncol = 1, 
      nrow = 2,
      padding = unit(0.3, "cm"),
      top = textGrob(
        "Damage Profiler Analysis (Part 1)", 
        gp = gpar(fontsize = 20, fontface = "bold", col = "#2c3e50"),
        just = "top"
      )
    )
    all_pages[[length(all_pages) + 1]] <- damage_page1
    
    # Second damage page (if more plots)
    if (length(damage_grobs) > 2) {
      damage_page2 <- arrangeGrob(
        grobs = damage_grobs[3:length(damage_grobs)], 
        ncol = 1,
        padding = unit(0.3, "cm"),
        top = textGrob(
          "Damage Profiler Analysis (Part 2)", 
          gp = gpar(fontsize = 20, fontface = "bold", col = "#2c3e50"),
          just = "top"
        )
      )
      all_pages[[length(all_pages) + 1]] <- damage_page2
    }
  } else {
    # Only 1 plot - use full page
    damage_page <- arrangeGrob(
      grobs = damage_grobs, 
      ncol = 1,
      padding = unit(0.3, "cm"),
      top = textGrob(
        "Damage Profiler Analysis", 
        gp = gpar(fontsize = 20, fontface = "bold", col = "#2c3e50"),
        just = "top"
      )
    )
    all_pages[[length(all_pages) + 1]] <- damage_page
  }
}

# ------------------------- EXPORT -------------------------
if (length(all_pages) == 0) {
  message("[WARNING] No plots available for ", sample, " - ", pathogen)
  quit(status = 0)
}

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

# Save as multi-page PDF with larger page size
pdf(output_path, width = 14, height = 10, onefile = TRUE, paper = "a4r")
for (i in seq_along(all_pages)) {
  page_obj <- all_pages[[i]]
  if (inherits(page_obj, "ggplot")) {
    # ggplot objects must be printed to render on the active device
    print(page_obj)
  } else {
    grid.draw(page_obj)
  }
  # Only add new page if not the last page
  if (i < length(all_pages)) {
    grid.newpage()
  }
}
dev.off()

message("[SUCCESS] Generated pathogen report: ", output_path)
