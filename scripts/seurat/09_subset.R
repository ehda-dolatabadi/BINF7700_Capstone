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
cores   <- args[4]

# Parameters
idents <- 14
dims <- 30
res <- 1.5
metric  <- "cosine"
seed    <- 777

significance <- 0.05
regulation <- 0.25
enrichment <- 0.1
top <- 10

# Set up parallelization
options(future.globals.maxSize = 16000 * 1024^2)  # 16 GB
plan("multicore", workers = as.numeric(cores))

# Load PCA object
obj <- readRDS(input)

# Subset to cluster
obj <- subset(obj, idents = idents)

# Switch back to RNA
DefaultAssay(obj) <- "RNA"

# Re-run SCTransform from scratch on the subset
obj <- SCTransform(
  obj,
  assay = "RNA",
  new.assay.name = "SCT",
)

# Set to SCT
DefaultAssay(obj) <- "SCT"

# Re-run PCA
obj <- RunPCA(obj, npcs = 50)

# Prepare SCT assay
obj <- PrepSCTFindMarkers(obj)

# Re-run clustering
obj <- FindNeighbors(obj, dims = 1:dims)
obj <- FindClusters(obj, res = res, algorithm = 4)

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_clustered.rds")))

# Re-run UMAP
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

# Find markers
markers <- FindAllMarkers(
  obj,
  only.pos = FALSE,
)

# Order markers
markers <- markers[order(markers$cluster, markers$p_val_adj, -markers$avg_log2FC), ]

# Save markers
outfile <- file.path(outdir, paste0(id, "_markers.tsv"))
write.table(markers, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)

# Filter
markers <- markers %>%
  mutate(direction = ifelse(avg_log2FC >= 0, "up", "down")) %>%
  filter(
    p_val_adj < significance,
    abs(avg_log2FC) >= regulation,
    abs(pct.1 - pct.2) >= enrichment
  )

# Rank and pick top rankings per cluster and direction
top <- markers %>%
  group_by(cluster, direction) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)), desc(abs(pct.1 - pct.2))) %>%
  slice_head(n = top) %>%
  ungroup()

# Save filtered markers
outfile <- file.path(outdir, paste0(id, "_markers_filtered.tsv"))
write.table(top, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
