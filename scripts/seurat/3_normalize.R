#!/usr/bin/env Rscript
# Purpose: SCTransform normalization on a filtered Seurat object
# Usage: Rscript 3_normalize.R <id> <input> <outdir>

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)

id	<- args[1]
input	<- args[2]
outdir	<- args[3]

# Load filtered Seurat object from QC step
obj <- readRDS(input)

# SCTransform
obj <- SCTransform(
  object = obj,
  assay = "RNA",
  new.assay.name = "SCT",
  verbose = TRUE
)

# Set default assay to SCT for downstream steps
DefaultAssay(obj) <- "SCT"

# Save outputs
saveRDS(obj, file = file.path(outdir, paste0(id, "_normalized.rds")))

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

print(sessionInfo())
