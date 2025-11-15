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
options(future.globals.maxSize = 16000 * 1024^2)  # 16 GB
plan("multicore", workers = as.numeric(cores))

# Load normalized objects
obj_list <- lapply(paths, readRDS)

# Integration features
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)

# Find anchors and integrate (batch correction)
anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  normalization.method = "SCT",
  anchor.features = features
)

obj <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

# Set default assay to integrated for downstream steps
DefaultAssay(obj) <- "integrated"

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_integrated.rds")))

# Summary
n_cells_total <- sum(sapply(obj_list, ncol))
write.table(
  data.frame(
    id				= id,
    n_samples			= length(obj_list),
    n_cells_total		= n_cells_total,
    n_features_integrated	= nrow(obj),
    default_assay		= DefaultAssay(obj),
    n_integration_features	= length(features)
  ),
  file = file.path(outdir, paste0(id, "_integrated.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
