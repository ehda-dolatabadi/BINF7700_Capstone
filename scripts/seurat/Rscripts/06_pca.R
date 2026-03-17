#!/usr/bin/env Rscript
# Script: 06_pca.R
# Purpose: Perform principal component analysis for dimensionality reduction
# Description: Runs PCA on integrated data and generates an elbow plot to help
#              determine the optimal number of PCs for downstream analysis
# Usage: Rscript 06_pca.R <id> <output> <input> <npcs>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(future)
})

set.seed(777)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]
npcs    <- as.numeric(args[4])

# Set up parallelization
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 56)

# Load integrated object
obj <- readRDS(input)

# PCA
obj <- RunPCA(obj, npcs = npcs)

# Clear scale.data (only needed for pca)
LayerData(obj, assay = DefaultAssay(obj), layer = "scale.data") <- NULL

# PCA info
pca_stdev <- obj@reductions$pca@stdev
total_var <- sum(pca_stdev^2)
pct_var_pc1 <- (pca_stdev[1]^2 / total_var) * 100

# Elbow plot
png(file.path(outdir, paste0(id, "_pca_elbow.png")), width=1600, height=1200, res=150)
print(ElbowPlot(obj, ndims = npcs) +
	geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
	labs(
		title = "PCA Elbow Plot",
		x = "Principal Component",
		y = "Standard Deviation"
	) +
	theme(
		plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
		axis.title = element_text(size = 14, face = "bold"),
		axis.text = element_text(size = 12),
		legend.title = element_text(size = 12, face = "bold"),
		legend.text = element_text(size = 11)
	))
dev.off()

# Save object
saveRDS(obj, file = file.path(dirname(outdir), paste0(id, "_processed.rds")))

# Summary table
write.table(
  data.frame(
    id,
    npcs,
    pct_var_pc1
  ),
  file = file.path(outdir, paste0(id, "_pca_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
