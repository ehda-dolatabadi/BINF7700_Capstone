#!/usr/bin/env Rscript
# Script: 04_normalize.R
# Purpose: Normalize filtered data using SCTransform
# Description: Applies SCTransform normalization to stabilize variance and identify
#              highly variable features for downstream analysis
# Usage: Rscript 04_normalize.R <id> <output> <input>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
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
plan("multicore", workers = 4)

# Load filtered Seurat object
obj <- readRDS(input)

# Apply SCTransform normalization
obj <- SCTransform(
  object = obj,
  assay = "RNA",
  new.assay.name = "SCT",
  ncells = 5000,		# default
  variable.features.n = 3000,	# default
  seed.use = seed
)
n_variable_features <- length(VariableFeatures(obj))

# Set SCT as default assay for downstream analysis
DefaultAssay(obj) <- "SCT"

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_normalized.rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells = ncol(obj),
    n_features = nrow(obj),
    n_variable_features
  ),
  file = file.path(outdir, paste0(id, "_normalization_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
