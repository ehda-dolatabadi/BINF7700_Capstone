#!/usr/bin/env Rscript
# Script: 07_cluster.R
# Purpose: Perform graph-based clustering using the Leiden algorithm
# Description: Builds nearest neighbor graph and identifies cell clusters using
#              Leiden community detection algorithm
# Usage: Rscript 07_cluster.R <id> <output> <input> <dims> <res>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(future)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(271)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]
dims    <- as.numeric(args[4])
res     <- as.numeric(args[5])

# Set up parallelization
options(future.seed = TRUE)
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 56)

# Load PCA object
obj <- readRDS(input)

# Neighbors and clustering (Leiden)
# FindNeighbors uses Annoy (C-level RNG) — set.seed() does not control it;
# switch to nn.method="rann" for R-level reproducibility
set.seed(271)
obj <- FindNeighbors(obj, dims = 1:dims)
obj <- FindClusters(
  obj,
  res = res,
  algorithm = 4,
  leiden_method = "igraph",
  random.seed = 271  # default random.seed: 0 (warns and resets to 1)
)

# Save object
saveRDS(obj, file = file.path(dirname(outdir), paste0(id, "_processed.rds")))

# Summary table
cluster_counts <- table(obj$seurat_clusters)
cluster_sizes <- paste(names(cluster_counts), cluster_counts, sep = ":", collapse = "; ")

write.table(
  data.frame(
    id,
    n_clusters = length(unique(obj$seurat_clusters)),
    cells_per_cluster = cluster_sizes,
    dims,
    res
  ),
  file = file.path(outdir, paste0(id, "_clustering_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
