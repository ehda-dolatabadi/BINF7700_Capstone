#!/usr/bin/env Rscript
# Purpose: Run UMAP on clustered Seurat object
# Usage: Rscript 7_umap.R <id> <output> <input>

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
dims	<- 10
metric	<- "cosine"
seed	<- 777

# Set up parallelization
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 64)

# Load object
obj <- readRDS(input)

# UMAP
obj <- RunUMAP(
  object       = obj,
  dims         = 1:dims,
  metric       = metric,
  seed.use     = seed
)

# Plot
png(file.path(outdir, paste0(id, "_07_umap.png")), width = 1600, height = 1200)
print(DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 0.5))
dev.off()

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_07_umap.rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells = ncol(obj),
    n_features = nrow(obj),
    n_clusters = length(unique(obj$seurat_clusters)),
    dims,
    metric,
    seed
  ),
  file = file.path(outdir, paste0(id, "_07_umap_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
