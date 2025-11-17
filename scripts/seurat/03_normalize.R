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
cores	<- args[4]

# Set up parallelization
options(future.globals.maxSize = 16000 * 1024^2)  # 16 GB
plan("multicore", workers = as.numeric(cores))

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

# Set default assay to SCT for downstream steps
DefaultAssay(obj) <- "SCT"

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_normalized.rds")))

# Summary
write.table(
  data.frame(
    id			= id,
    n_cells		= ncol(obj),
    n_features		= nrow(obj),
    default_assay	= DefaultAssay(obj),
    n_var_features	= length(VariableFeatures(obj))
  ),
  file = file.path(outdir, paste0(id, "_normalized.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
