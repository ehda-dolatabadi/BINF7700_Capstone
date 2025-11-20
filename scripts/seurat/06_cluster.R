#!/usr/bin/env Rscript
# Purpose: Run graph-based clustering on PCA Seurat object
# Usage: Rscript 6_cluster.R <id> <output> <input>

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

# Parameters
dims <- 10
res <- 0.5

# Set up parallelization
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 64)

# Load PCA object
obj <- readRDS(input)

# Neighbors and clustering (Leiden)
obj <- FindNeighbors(obj, dims = 1:dims)
obj <- FindClusters(obj, res = res, algorithm = 4)

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_06_clustered.rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells = ncol(obj),
    n_features = nrow(obj),
    dims,
    res,
    n_clusters = length(unique(obj$seurat_clusters))
  ),
  file = file.path(outdir, paste0(id, "_06_clustering_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
