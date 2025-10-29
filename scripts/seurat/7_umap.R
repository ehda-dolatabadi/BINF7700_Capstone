#!/usr/bin/env Rscript
# Purpose: Run UMAP on clustered Seurat object
# Usage: Rscript 7_umap.R <id> <output> <input>

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]

# Parameters
dims	<- 20
metric	<- "cosine"
seed	<- 777

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
png(file.path(outdir, paste0(id, "_umap.png")), width = 1600, height = 1200)
print(DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 0.5))
dev.off()

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_umap.rds")))

# Summary
write.table(
  data.frame(
    id               = id,
    n_cells          = ncol(obj),
    n_features       = nrow(obj),
    dims_used        = dims,
    umap_metric      = metric,
    umap_seed        = seed
  ),
  file = file.path(outdir, paste0(id, "_umap.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
