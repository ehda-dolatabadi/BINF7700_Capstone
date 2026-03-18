#!/usr/bin/env Rscript
# Script: 05_reintegrate_subset.R
# Purpose: Re-normalize and re-integrate a subset of cells
# Description: Splits subset by sample, re-runs SCTransform on each, and performs
#              integration (skips integration if any sample has <30 cells)
# Usage: Rscript 05_reintegrate_subset.R <id> <outdir> <input>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
})

seed <- 777
set.seed(seed)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]

# Configure parallel processing
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 7)

# Load subset Seurat object
obj <- readRDS(input)
n_cells_before <- ncol(obj)
n_features_before <- nrow(obj)

# Switch to RNA assay for re-normalization
DefaultAssay(obj) <- "RNA"

# clear extra assays
# obj[["SCT"]] <- NULL
# obj[["integrated"]] <- NULL

# Split by sample and check if integration is feasible
obj_list <- SplitObject(obj, split.by = "orig.ident")
cell_counts <- sapply(obj_list, ncol)
samples <- names(cell_counts)[cell_counts < 30]

if (length(samples) > 0) {
  cat("Some samples had less than 30 cells:", paste(samples, collapse = ", "), "\n")
  cat("Running SCTransform on whole object without re-integration\n")
  
  # Run SCTransform on the whole object (not split)
  obj <- SCTransform(
    obj,
    assay = "RNA",
    new.assay.name = "SCT",
  )
  
  # Set default to SCT (not integrated, since integration was skipped)
  DefaultAssay(obj) <- "SCT"
  
  # Save object
  saveRDS(obj, file = file.path(dirname(outdir), paste0(id, "_integrated.rds")))

  # Summary table
  write.table(
    data.frame(
      id,
      n_cells_before,
      n_cells_after = ncol(obj),
      n_features_before,
      n_features_after = nrow(obj),
      integration_performed = FALSE
    ),
    file = file.path(outdir, paste0(id, "_integration_summary.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  quit(save = "no", status = 0)
}

# Re-run SCTransform on each sample separately
obj_list <- lapply(obj_list, function(x) {
  SCTransform(x, assay = "RNA", new.assay.name = "SCT")
})

# Set to SCT assay for each object in the list
obj_list <- lapply(obj_list, function(x) {
  DefaultAssay(x) <- "SCT"
  x
})

# Integration features
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)

# Find anchors and integrate (batch correction)
anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  normalization.method = "SCT",
  anchor.features = features
)
n_anchors <- nrow(anchors@anchors)

obj <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

# Set default assay to integrated for downstream steps
DefaultAssay(obj) <- "integrated"

# Clear scale.data (only needed for integrate)
# LayerData(obj, assay = "SCT", layer = "scale.data") <- NULL

# Correct timepoints order
timepoint_order <- c("control", "3h", "24h", "72h", "7dpa", "14dpa", "22dpa", "33dpa")
obj$orig.ident <- factor(obj$orig.ident, levels = timepoint_order)

# Save object
saveRDS(obj, file = file.path(dirname(outdir), paste0(id, "_integrated.rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells_before,
    n_cells_after = ncol(obj),
    n_features_before,
    n_features_after = nrow(obj)
  ),
  file = file.path(outdir, paste0(id, "_integration_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
