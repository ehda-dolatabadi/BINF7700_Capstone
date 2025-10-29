#!/usr/bin/env Rscript
# Purpose: Generate UMAP plot for clustered Seurat object
# Usage: Rscript 8_plot_umap.R <id> <output> <input>

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]

# Load object
obj <- readRDS(input)

# Plot
png(file.path(outdir, paste0(id, "_umap.png")), width = 1600, height = 1200)
print(DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 0.5))
dev.off()

# Summary
write.table(
  data.frame(
    id		= id,
    n_cells	= ncol(obj),
    n_features	= nrow(obj),
    n_clusters	= length(levels(obj$seurat_clusters))
  ),
  file = file.path(outdir, paste0(id, "_umap_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(sessionInfo())
