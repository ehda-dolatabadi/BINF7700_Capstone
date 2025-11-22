#!/usr/bin/env Rscript
# Purpose: SCTransform normalization on a filtered Seurat object
# Usage: Rscript 3_normalize.R <id> <output> <input>

suppressPackageStartupMessages({
  library(Seurat)
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

# Load filtered object
obj <- readRDS(input)

# SCTransform
obj <- SCTransform(
  object = obj,
  assay = "RNA",
  new.assay.name = "SCT",
  ncells = 5000,		# default
  variable.features.n = 3000	# default
)
n_variable_features <- length(VariableFeatures(obj))

# Set default assay to SCT for downstream steps
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
