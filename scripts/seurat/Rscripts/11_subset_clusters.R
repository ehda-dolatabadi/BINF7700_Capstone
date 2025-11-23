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
idents	<- args[4]

# Set up parallelization
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 8)

# Load clustered object
obj <- readRDS(input)
n_cells_before <- ncol(obj)
n_features_before <- nrow(obj)

# Subset to cluster
obj <- subset(obj, idents = idents)

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_00_subset.rds")))

# Summary table
write.table(
  data.frame(
    id,
    idents,
    n_cells_before,
    n_cells_after = ncol(obj),
    n_features_before,
    n_features_after = nrow(obj)
  ),
  file = file.path(outdir, paste0(id, "_00_subsetting_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
