#!/usr/bin/env Rscript
# Purpose: Integrate multiple SCT-normalized Seurat objects
# Usage: Rscript 4_integrate.R <id> <outdir> <sample1.rds> <sample2.rds> ...

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id	<- args[1]
outdir	<- args[2]
paths	<- args[3:length(args)]

# Set up parallelization
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 8)

# Load normalized objects
obj_list <- lapply(paths, readRDS)
n_samples <- length(obj_list)
n_features_per_sample <- sapply(obj_list, nrow)
n_features_min <- min(n_features_per_sample)
n_features_max <- max(n_features_per_sample)
n_features_mean <- mean(n_features_per_sample)

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

# Correct timepoints order
timepoint_order <- c("control", "3h", "24h", "72h", "7dpa", "14dpa", "22dpa", "33dpa")
obj$orig.ident <- factor(obj$orig.ident, levels = timepoint_order)

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_integrated.rds")))

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
