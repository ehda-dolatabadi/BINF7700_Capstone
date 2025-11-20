#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]

# Set up parallelization
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 8)

# Load subset object
obj <- readRDS(input)
n_cells_before <- ncol(obj)
n_features_before <- nrow(obj)

# Switch back to RNA assay
DefaultAssay(obj) <- "RNA"

# Re-run SCTransform on the subset
obj <- SCTransform(
  obj,
  assay = "RNA",
  new.assay.name = "SCT",
)

# Set to SCT assay
DefaultAssay(obj) <- "SCT"

# Split by sample
obj_list <- SplitObject(obj, split.by = "orig.ident")

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
saveRDS(obj, file = file.path(outdir, paste0(id, "_04_integrated.rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells_before,
    n_cells_after = ncol(obj),
    n_features_before,
    n_features_after = nrow(obj)
  ),
  file = file.path(outdir, paste0(id, "_04_integration_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
