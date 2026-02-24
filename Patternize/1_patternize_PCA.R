#!/usr/bin/env Rscript
###############################################################################
# 01_patternize_PCA.R
#
# Create PCA scores from wing images using patternize:
#   - Load images from a folder
#   - (Optional) filter to a subset list
#   - Choose a reference image
#   - Choose an RGB seed (interactive sampling or fixed RGB)
#   - Register patterns with patRegRGB()
#   - Run PCA with patPCA()
#   - Write PCA scores to CSV
###############################################################################

suppressPackageStartupMessages({
  library(patternize)
  library(tools)
  library(raster)
  library(viridis)
  library(jpeg)
  library(png)
  library(factoextra)
})

# --- USER SETTINGS ------------------------------------------------------------
img_dir <- "../2_colour_corrected_cropped/HW/"  # folder containing images
img_ext <- ".jpeg"                              # ".jpeg", ".jpg", ".png", etc.

# Optional: subset list
#   - NULL (use all), OR a vector of filenames, OR a text file with one filename/line
subset <- NULL
# subset <- "specimens_to_use.txt"

# Reference specimen ID = filename without extension (set NULL to use first image)
ref_id <- NULL

# RGB seed selection
use_interactive_rgb <- TRUE
fixed_rgb <- c(114, 96, 58)  # used only if use_interactive_rgb = FALSE

# patRegRGB parameters
resample_factor <- 1
col_offset      <- 0.05
remove_bg       <- 100
plot_registration <- "stack"   # "stack" or FALSE

# PCA settings
pcs_to_export <- 5
plot_pca      <- TRUE

# Output
out_pca_csv <- "01_PCA_scores.csv"

# --- Helpers -----------------------------------------------------------------
read_subset <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 1 && file.exists(x)) {
    vals <- trimws(readLines(x, warn = FALSE))
    vals <- vals[vals != ""]
    return(vals)
  }
  as.character(x)
}

# --- Load images --------------------------------------------------------------
if (!dir.exists(img_dir)) stop("img_dir does not exist: ", img_dir)

ext_pattern <- paste0("\\", img_ext, "$")
img_files <- list.files(img_dir, pattern = ext_pattern, full.names = FALSE)
if (length(img_files) == 0) stop("No files with extension ", img_ext, " found in ", img_dir)

subset_files <- read_subset(subset)
if (!is.null(subset_files)) {
  img_files <- img_files[img_files %in% subset_files]
  if (length(img_files) == 0) stop("After filtering, no files remain. Check subset names/extensions.")
}

id_list <- file_path_sans_ext(img_files)
image_list <- makeList(id_list, type = "image", prepath = img_dir, extension = img_ext)

if (is.null(ref_id)) {
  ref_id <- names(image_list)[1]
  message("ref_id not provided; using first image as reference: ", ref_id)
}
if (!ref_id %in% names(image_list)) stop("ref_id not found in image_list: ", ref_id)

# --- RGB seed -----------------------------------------------------------------
if (use_interactive_rgb) {
  message("Sampling RGB from reference image: ", ref_id)
  rgb_seed <- sampleRGB(image_list[[ref_id]], resampleFactor = NULL)
} else {
  if (length(fixed_rgb) != 3) stop("fixed_rgb must be length 3.")
  rgb_seed <- fixed_rgb
}
message("Using RGB seed: ", paste(rgb_seed, collapse = ", "))

# --- Register patterns --------------------------------------------------------
raster_list_regRGB <- patRegRGB(
  image_list,
  image_list[[ref_id]],
  rgb_seed,
  resampleFactor = resample_factor,
  colOffset      = col_offset,
  removebg       = remove_bg,
  plot           = plot_registration
)

# --- PCA ----------------------------------------------------------------------
pca_out <- patPCA(
  raster_list_regRGB,
  popList     = list(names(image_list)),
  colList     = "red",
  plot        = plot_pca,
  plotType    = "points",
  plotChanges = TRUE,
  PCx         = 1,
  PCy         = 2,
  main        = "PCA of Wing Patterns",
  normalized  = TRUE
)

scores <- pca_out[[3]]$x
specimen_ids <- pca_out[[2]]$sampleID

n_pc <- min(pcs_to_export, ncol(scores))
pc_df <- data.frame(SpecimenID = specimen_ids, scores[, seq_len(n_pc), drop = FALSE])
colnames(pc_df)[-1] <- paste0("PC", seq_len(n_pc))

write.csv(pc_df, out_pca_csv, row.names = FALSE)
message("Done. Wrote PCA scores to: ", out_pca_csv)
