#!/usr/bin/env Rscript
# Generate a comprehensive per-sample PDF report

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(readr)
  library(readxl)
})

# NOTE: Animal icons previously used the magick package to embed SVGs into the
# PDF header. On some systems, when many sample-report jobs run in parallel,
# this triggered low-level "free(): double free detected in tcache 2" errors
# from underlying C libraries. To keep the pipeline robust, we disable icons
# entirely here (text-only header).
USE_ICONS <- FALSE

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript generate_sample_report.R <sample> <host_mtdna_summary> <pathogen_summary> <output_pdf> [animal_icons_dir]")
}

sample <- args[1]
host_mtdna_file <- args[2]
pathogen_summary_file <- args[3]
output_pdf <- args[4]
animal_icons_dir <- if (length(args) >= 5) args[5] else "config/animal_icons"  # Default directory

# Read summary files
# Note: 'sample' argument is the biological sample (bio_sample), not PCR-level
host_mtdna <- tryCatch({
  # FIXED: read_excel() doesn't have show_col_types parameter (that's for read_csv())
  df <- read_excel(host_mtdna_file, sheet = "Sample_level")
  # Normalize column names: "bio_sample" -> "Sample" for consistency
  # The host/mtDNA summary uses "bio_sample" as the column name
  if ("bio_sample" %in% colnames(df)) {
    if (!"Sample" %in% colnames(df)) {
      df$Sample <- df$bio_sample
    }
  }
  # Also normalize "species" to "Species" if needed
  if ("species" %in% colnames(df) && !"Species" %in% colnames(df)) {
    df$Species <- df$species
  }
  df
}, error = function(e) {
  message("[WARNING] Could not read host/mtDNA summary: ", e$message)
  NULL
})

pathogen_summary <- tryCatch({
  # Detect file type: Excel vs CSV
  if (grepl("\\.xlsx?$", pathogen_summary_file, ignore.case = TRUE)) {
    read_excel(pathogen_summary_file)
  } else {
    read_csv(pathogen_summary_file, show_col_types = FALSE)
  }
}, error = function(e) {
  message("[WARNING] Could not read pathogen summary: ", e$message)
  NULL
})

# Filter to current biological sample
# The 'sample' argument is the biological sample name (bio_sample)
if (!is.null(host_mtdna) && nrow(host_mtdna) > 0) {
  # After normalization, we should have "Sample" column (from "bio_sample")
  # But also check for "bio_sample" in case normalization didn't happen
  sample_col <- if ("Sample" %in% colnames(host_mtdna)) "Sample" else if ("bio_sample" %in% colnames(host_mtdna)) "bio_sample" else NULL
  if (!is.null(sample_col)) {
    # Store available samples before filtering for error message
    available_samples <- unique(host_mtdna[[sample_col]])
    # Filter by biological sample name
    host_mtdna <- host_mtdna[host_mtdna[[sample_col]] == sample, ]
    if (nrow(host_mtdna) == 0) {
      message("[WARNING] Biological sample '", sample, "' not found in host/mtDNA summary.")
      message("[INFO] Available biological samples: ", paste(head(available_samples, 10), collapse = ", "))
      host_mtdna <- NULL
    }
  } else {
    message("[WARNING] Could not find sample column in host/mtDNA summary. Available columns: ", 
            paste(colnames(host_mtdna), collapse = ", "))
    host_mtdna <- NULL
  }
} else if (!is.null(host_mtdna)) {
  message("[WARNING] host/mtDNA summary is empty")
  host_mtdna <- NULL
}

if (!is.null(pathogen_summary) && nrow(pathogen_summary) > 0) {
  if ("Sample" %in% colnames(pathogen_summary)) {
    pathogen_summary <- pathogen_summary[pathogen_summary$Sample == sample, ]
    if (nrow(pathogen_summary) == 0) {
      message("[INFO] No pathogens detected for sample '", sample, "'")
      pathogen_summary <- NULL
    }
  } else {
    message("[WARNING] Could not find 'Sample' column in pathogen summary. Available columns: ", 
            paste(colnames(pathogen_summary), collapse = ", "))
    pathogen_summary <- NULL
  }
} else if (!is.null(pathogen_summary)) {
  message("[WARNING] Pathogen summary is empty")
  pathogen_summary <- NULL
}

# Normalize pathogen coverage column: create a unified 'Coverage' column if needed
if (!is.null(pathogen_summary) && nrow(pathogen_summary) > 0) {
  if (!"Coverage" %in% colnames(pathogen_summary)) {
    # Look for any column that looks like pathogen coverage
    cov_candidates <- colnames(pathogen_summary)[grepl("cover", tolower(colnames(pathogen_summary)))]
    if (length(cov_candidates) > 0) {
      # Use the first candidate as Coverage
      pathogen_summary$Coverage <- pathogen_summary[[cov_candidates[1]]]
    }
  }
}


# Function to find and load animal icon SVG
load_animal_icon <- function(species_name, icons_dir) {
  # Early exit if icons are disabled or magick is unavailable
  if (!USE_ICONS) {
    return(NULL)
  }

  if (is.null(species_name) || is.na(species_name) || species_name == "NA" || !dir.exists(icons_dir)) {
    return(NULL)
  }
  
  # Common name mapping for scientific names (you can extend this)
  common_name_map <- list(
    "ovis aries" = "sheep",
    "bos taurus" = "cattle",
    "bos indicus" = "cattle",
    "equus caballus" = "horse",
    "sus scrofa" = "pig",
    "capra hircus" = "goat",
    "gallus gallus" = "chicken",
    "canis familiaris" = "dog",
    "felis catus" = "cat",
    "homo sapiens" = "human"
  )
  
  # Normalize species name for matching: lowercase, replace spaces with underscores or hyphens
  normalize_name <- function(name) {
    name <- tolower(trimws(name))
    # Try multiple variations
    variations <- c(
      name,
      gsub(" ", "_", name),
      gsub(" ", "-", name),
      gsub("_", "-", name),
      gsub("-", "_", name)
    )
    return(unique(variations))
  }
  
  # Get both scientific and common name variations
  species_variations <- normalize_name(species_name)
  
  # Add common name if mapping exists
  species_lower <- tolower(trimws(species_name))
  if (species_lower %in% names(common_name_map)) {
    common_name <- common_name_map[[species_lower]]
    species_variations <- c(species_variations, normalize_name(common_name))
  }
  
  # Look for SVG files matching any variation
  svg_files <- list.files(icons_dir, pattern = "\\.svg$", ignore.case = TRUE, full.names = TRUE)
  
  for (variation in species_variations) {
    # Try exact match first
    matches <- grep(paste0("^", variation, "\\.svg$"), basename(svg_files), ignore.case = TRUE)
    if (length(matches) > 0) {
      svg_path <- svg_files[matches[1]]
      tryCatch({
        # Read SVG and convert to raster image (magick can read SVG directly)
        img <- image_read(svg_path)
        # Resize to reasonable size for PDF
        img <- image_resize(img, "200x200")
        return(img)
      }, error = function(e) {
        message("[WARNING] Could not load SVG icon: ", svg_path, " - ", e$message)
        return(NULL)
      })
    }
    
    # Try partial match (species name contained in filename)
    matches <- grep(variation, basename(svg_files), ignore.case = TRUE)
    if (length(matches) > 0) {
      svg_path <- svg_files[matches[1]]
      tryCatch({
        img <- image_read(svg_path)
        img <- image_resize(img, "200x200")
        return(img)
      }, error = function(e) {
        message("[WARNING] Could not load SVG icon: ", svg_path, " - ", e$message)
        return(NULL)
      })
    }
  }
  
  return(NULL)
}

# Get species name for icon lookup
species_name <- NULL
if (!is.null(host_mtdna) && nrow(host_mtdna) > 0) {
  species_name <- host_mtdna$Species[1]
  if (is.na(species_name) || species_name == "NA") {
    # Try lowercase "species" column
    if ("species" %in% colnames(host_mtdna)) {
      species_name <- host_mtdna$species[1]
    }
  }
}

# Load animal icon if available
animal_icon <- NULL
if (!is.null(species_name)) {
  animal_icon <- load_animal_icon(species_name, animal_icons_dir)
}

# Create PDF - single page only (landscape A4)
pdf(output_pdf, width = 11.69, height = 8.27, onefile = TRUE, paper = "a4r")

# Single combined page: Sample overview + pathogen summary
if (!is.null(host_mtdna) && nrow(host_mtdna) > 0) {
  # Helper function to safely extract and format values (case-insensitive, with a bit of flexibility)
  get_val <- function(col_name, default = "NA", fmt = NULL) {
    # Try exact match first, then case-insensitive, then a few simple heuristics
    matching_col <- NULL
    col_lower <- tolower(col_name)
    for (col in colnames(host_mtdna)) {
      if (col == col_name || tolower(col) == col_lower) {
        matching_col <- col
        break
      }
    }
    
    # Extra handling for mtDNA coverage if exact name not found
    if (is.null(matching_col) && col_lower == "mtdna_coverage") {
      for (col in colnames(host_mtdna)) {
        col_l <- tolower(col)
        # Accept a broad set of mtDNA coverage column names:
        # must contain "cov" and either "mtdna", "mt", or "mito"
        if (grepl("cov", col_l) && (grepl("mtdna", col_l) || grepl("\\bmt", col_l) || grepl("mito", col_l))) {
          matching_col <- col
          break
        }
      }
    }
    
    if (is.null(matching_col)) return(default)
    val <- host_mtdna[[matching_col]][1]
    if (is.na(val) || val == "NA" || val == "" || is.null(val)) return(default)
    if (!is.null(fmt)) {
      tryCatch({
        num_val <- as.numeric(val)
        if (is.na(num_val)) return(default)
        return(fmt(num_val))
      }, error = function(e) return(default))
    }
    return(as.character(val))
  }
  
  overview_data <- data.frame(
    Metric = c("Species", "Endogenous %", "Host Coverage", "mtDNA Coverage", 
               "Total Reads", "Q30 Reads", "Duplication Rate %"),
    Value = c(
      get_val("species", "NA"),  # Note: lowercase "species" in the file
      get_val("host_endogenous_pct", "NA", function(x) paste0(round(x, 2), "%")),
      get_val("host_coverage", "NA", function(x) paste0(round(x, 4), "X")),
      get_val("mtdna_coverage", "NA", function(x) paste0(round(x, 4), "X")),
      get_val("raw_reads", "NA", function(x) format(x, big.mark = ",")),
      get_val("q30_reads", "NA", function(x) format(x, big.mark = ",")),
      get_val("duplication_rate", "NA", function(x) paste0(round(x, 2), "%"))
    )
  )
  
  overview_table <- tableGrob(
    overview_data, 
    rows = NULL,
    theme = ttheme_minimal(
      core = list(fg_params = list(fontsize = 9)),
      colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
    )
  )
}

# Page 2: Combined single-page layout (overview + pathogen summary)

# Build pathogen summary table (if any pathogens for this sample)
pathogen_table <- NULL
if (!is.null(pathogen_summary) && nrow(pathogen_summary) > 0) {
  available_cols <- c(
    "Pathogen", "Score", "Coverage", "Mapped_reads",
    "ANI", "Relative_entropy", "Breadth_ratio",
    "Damage_5p_CtoT", "Damage_3p_GtoA"
  )
  pathogen_table_cols <- available_cols[available_cols %in% colnames(pathogen_summary)]
  if (length(pathogen_table_cols) > 0) {
    pathogen_table_data <- pathogen_summary[, pathogen_table_cols]
    
    # Rename columns for display
    display_names <- c(
      "Pathogen" = "Pathogen",
      "Score" = "Score",
      "Coverage" = "Coverage",
      "Mapped_reads" = "Mapped Reads",
      "ANI" = "ANI (%)",
      "Relative_entropy" = "Entropy",
      "Breadth_ratio" = "Breadth Ratio",
      "Damage_5p_CtoT" = "5' C→T (%)",
      "Damage_3p_GtoA" = "3' G→A (%)"
    )
    colnames(pathogen_table_data) <- display_names[colnames(pathogen_table_data)]
    
    # Format numeric columns
    if ("Coverage" %in% colnames(pathogen_table_data)) {
      pathogen_table_data$Coverage <- sapply(pathogen_table_data$Coverage, function(x) {
        if (is.na(x) || x == "NA") return("NA")
        # Strip non-numeric characters (e.g., "X", spaces) before conversion
        x_num <- suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", as.character(x))))
        if (is.na(x_num)) return("NA")
        paste0(round(x_num, 4), "X")
      })
    }
    
    if ("ANI (%)" %in% colnames(pathogen_table_data)) {
      pathogen_table_data$`ANI (%)` <- sapply(pathogen_table_data$`ANI (%)`, function(x) {
        if (is.na(x) || x == "NA") return("NA")
        # Strip any non-numeric characters (e.g. "99.15%") before conversion
        x_num <- suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", as.character(x))))
        if (is.na(x_num)) return("NA")
        paste0(round(x_num, 2), "%")
      })
    }
    
    pathogen_table <- tableGrob(
      pathogen_table_data,
      rows = NULL,
      theme = ttheme_minimal(
        core = list(fg_params = list(fontsize = 8)),
        colhead = list(fg_params = list(fontsize = 9, fontface = "bold"))
      )
    )
  }
}

# Create a compact single-page layout using two-column design
# Header row
header_grobs <- list()
if (!is.null(animal_icon)) {
  tryCatch({
    animal_icon_resized <- image_resize(animal_icon, "80x80")
    img_raster <- as.raster(animal_icon_resized)
    icon_grob <- rasterGrob(
      img_raster,
      width = unit(0.8, "inches"),
      height = unit(0.8, "inches"),
      interpolate = TRUE
    )
    header_grobs <- list(
      icon_grob,
      textGrob(
        paste0("Sample Report: ", sample),
        gp = gpar(fontsize = 18, fontface = "bold", col = "#2c3e50"),
        x = 0, just = "left"
      )
    )
  }, error = function(e) {
    header_grobs <<- list(
      textGrob(
        paste0("Sample Report: ", sample),
        gp = gpar(fontsize = 18, fontface = "bold", col = "#2c3e50")
      )
    )
  })
} else {
  header_grobs <- list(
    textGrob(
      paste0("Sample Report: ", sample),
      gp = gpar(fontsize = 18, fontface = "bold", col = "#2c3e50")
    )
  )
}

if (length(header_grobs) > 1) {
  header_row <- grid.arrange(
    grobs = header_grobs,
    ncol = 2,
    widths = c(0.08, 0.92)
  )
} else {
  header_row <- header_grobs[[1]]
}

subtitle_row <- textGrob(
  "PIGSTI Pipeline Summary",
  gp = gpar(fontsize = 11, col = "#7f8c8d")
)

# Main content: two columns side by side
left_column <- NULL
right_column <- NULL

# Left column: Sample Overview
if (!is.null(host_mtdna) && nrow(host_mtdna) > 0) {
  left_column <- grid.arrange(
    grobs = list(
      textGrob(
        "Sample Overview",
        gp = gpar(fontsize = 12, fontface = "bold", col = "#2c3e50")
      ),
      overview_table
    ),
    ncol = 1,
    heights = c(0.12, 0.88)
  )
}

# Right column: Pathogen Detection Summary
if (!is.null(pathogen_table)) {
  right_column <- grid.arrange(
    grobs = list(
      textGrob(
        "Pathogen Detection Summary",
        gp = gpar(fontsize = 12, fontface = "bold", col = "#2c3e50")
      ),
      pathogen_table
    ),
    ncol = 1,
    heights = c(0.12, 0.88)
  )
}

# Combine everything into single page
main_content <- NULL
if (!is.null(left_column) && !is.null(right_column)) {
  # Both columns available - side by side
  main_content <- grid.arrange(
    grobs = list(left_column, right_column),
    ncol = 2,
    widths = c(0.45, 0.55)
  )
} else if (!is.null(left_column)) {
  # Only left column
  main_content <- left_column
} else if (!is.null(right_column)) {
  # Only right column
  main_content <- right_column
}

# Final layout: header, subtitle, main content
if (!is.null(main_content)) {
  grid.arrange(
    grobs = list(
      header_row,
      subtitle_row,
      main_content
    ),
    ncol = 1,
    heights = c(0.08, 0.03, 0.89),
    padding = unit(0.5, "cm")
  )
} else {
  # Fallback if no content
  grid.arrange(
    grobs = list(
      header_row,
      subtitle_row
    ),
    ncol = 1,
    heights = c(0.5, 0.5),
    padding = unit(0.5, "cm")
  )
}

dev.off()

message("[SUCCESS] Generated sample report: ", output_pdf)

