#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(SingleCellExperiment)
  library(scDblFinder) 
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]

# Load mapped object
obj <- readRDS(input)

# Store initial counts
n_cells_before <- ncol(obj)
n_features_before <- nrow(obj)

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(obj)

# Run scDblFinder
sce <- scDblFinder(sce)

# Create vector of singlets
singlet_barcodes <- rownames(colData(sce))[sce$scDblFinder.class == "singlet"]

# Subset out doublets in the Seurat object
obj <- subset(obj, cells = singlet_barcodes)

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_DB_removed.rds")))

# Summary table
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
