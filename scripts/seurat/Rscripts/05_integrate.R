#!/usr/bin/env Rscript
# Script: 05_integrate.R
# Purpose: Integrate multiple normalized samples using Seurat's integration method
# Description: Performs batch correction across samples by finding integration anchors
#              and integrating data to enable cross-sample comparisons
# Usage: Rscript 05_integrate.R <id> <outdir> <sample1.rds> <sample2.rds> ...

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(future)
})

seed <- 777
set.seed(seed)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id	<- args[1]
outdir	<- args[2]
paths	<- args[3:length(args)]

# Configure parallel processing
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 7)

# Load all normalized Seurat objects
obj_list <- lapply(paths, readRDS)
n_samples <- length(obj_list)
n_features_per_sample <- sapply(obj_list, nrow)
n_features_min <- min(n_features_per_sample)
n_features_max <- max(n_features_per_sample)
n_features_mean <- mean(n_features_per_sample)

# Select integration features and prepare objects
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)

# Find integration anchors across samples
anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  normalization.method = "SCT",
  anchor.features = features
)
n_anchors <- nrow(anchors@anchors)

# Integrate data using anchors for batch correction
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
    n_samples,
    n_cells_total = ncol(obj),
    n_features_min,
    n_features_max,
    n_features_mean,
    n_features_integrated = nrow(obj),
    n_anchors
  ),
  file = file.path(outdir, paste0(id, "_integration_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
