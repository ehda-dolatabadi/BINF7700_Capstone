#!/usr/bin/env Rscript
# Script: 13_subset_clusters.R
# Purpose: Extract cells from specific cluster identities
# Description: Subsets the Seurat object to retain only cells from specified clusters
#              for focused downstream analysis
# Usage: Rscript 13_subset_clusters.R <id> <outdir> <input> <idents>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
})

set.seed(777)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]
idents	<- args[4]

# Set up parallelization
options(future.seed = TRUE)
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 56)

# Load clustered object
obj <- readRDS(input)
n_cells_before <- ncol(obj)
n_features_before <- nrow(obj[["SCT"]])

# Subset to cluster
obj <- subset(obj, idents = idents)

# Save object
saveRDS(obj, file = file.path(dirname(outdir), paste0(id, "_integrated.rds")))

# Summary table
write.table(
  data.frame(
    id,
    idents,
    n_cells_before,
    n_cells_after = ncol(obj),
    n_features_before,
    n_features_after = nrow(obj[["SCT"]])
  ),
  file = file.path(outdir, paste0(id, "_subsetting_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
