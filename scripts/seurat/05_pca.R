#!/usr/bin/env Rscript
# Purpose: Run PCA on an integrated Seurat object
# Usage: Rscript 5_pca.R <id> <output> <input>

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(future)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]
cores   <- args[4]

# Set up parallelization
options(future.globals.maxSize = 16000 * 1024^2)  # 16 GB
plan("multicore", workers = as.numeric(cores))

# Parameters
npcs <- 50

# Load integrated object
obj <- readRDS(input)

# PCA
obj <- RunPCA(obj, npcs = npcs)

# Elbow plot
png(file.path(outdir, paste0(id, "_elbow.png")), width=1200, height=720)
print(ElbowPlot(obj, ndims = npcs))
dev.off()

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_pca.rds")))

# Summary
stdev <- Stdev(obj, reduction = "pca")
pct_var <- (stdev^2) / sum(stdev^2) * 100

write.table(
  data.frame(
    id		= id,
    n_cells	= ncol(obj),
    n_features	= nrow(obj),
    n_pcs	= length(stdev),
    var_top10	= round(sum(pct_var[1:10]), 2),
    var_top30	= round(sum(pct_var[1:30]), 2),
    var_total	= round(sum(pct_var), 2)
  ),
  file = file.path(outdir, paste0(id, "_pca.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
