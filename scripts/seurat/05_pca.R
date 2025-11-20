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

# Set up parallelization
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 64)

# Parameters
npcs <- 50

# Load integrated object
obj <- readRDS(input)

# PCA
obj <- RunPCA(obj, npcs = npcs)

# PCA info
pca_stdev <- obj@reductions$pca@stdev
total_var <- sum(pca_stdev^2)
pct_var_pc1 <- (pca_stdev[1]^2 / total_var) * 100

# Elbow plot
png(file.path(outdir, paste0(id, "_05_pca_elbow.png")), width=1200, height=720)
print(ElbowPlot(obj, ndims = npcs) + 
	geom_vline(xintercept = 10, linetype = "dashed", color = "red"))
dev.off()

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_05_pca.rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells = ncol(obj),
    n_features = nrow(obj),
    npcs,
    pct_var_pc1
  ),
  file = file.path(outdir, paste0(id, "_05_pca_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
