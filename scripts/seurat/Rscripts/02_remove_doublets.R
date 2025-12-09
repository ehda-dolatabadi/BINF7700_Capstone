#!/usr/bin/env Rscript
# Script: 02_remove_doublets.R
# Purpose: Detect and remove doublets using scDblFinder
# Description: Converts Seurat object to SingleCellExperiment, runs doublet detection,
#              and filters out predicted doublets to retain only singlet cells
# Usage: Rscript 02_remove_doublets.R <id> <outdir> <input>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(SingleCellExperiment)
  library(scDblFinder)
})

set.seed(777)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]

# Load Seurat object
obj <- readRDS(input)

# Store initial counts
n_cells_before <- ncol(obj)
n_features_before <- nrow(obj)

# Convert to SingleCellExperiment format for scDblFinder
sce <- as.SingleCellExperiment(obj)

# Run doublet detection
sce <- scDblFinder(sce)

# Extract singlet cell barcodes
singlet_barcodes <- rownames(colData(sce))[sce$scDblFinder.class == "singlet"]

# Filter Seurat object to retain only singlets
obj <- subset(obj, cells = singlet_barcodes)

# Save filtered object
saveRDS(obj, file = file.path(outdir, paste0(id, "_DB_removed.rds")))

# Write summary statistics to file
write.table(
  data.frame(
    id,
    n_cells_before,
    n_cells_after = ncol(obj),
    n_features_before,
    n_features_after = nrow(obj)
  ),
  file = file.path(outdir, paste0(id, "_DB_removed_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
