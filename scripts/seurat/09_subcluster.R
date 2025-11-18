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
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 64)

# Load PCA object
obj <- readRDS(input)
n_cells_before <- ncol(obj)
n_features_before <- nrow(obj)

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
filtered <- markers %>%
  mutate(direction = ifelse(avg_log2FC >= 0, "up", "down")) %>%
  filter(
    p_val_adj < significance,
    abs(avg_log2FC) >= regulation,
    abs(pct.1 - pct.2) >= enrichment
  )

# Rank and pick top rankings per cluster and direction
filtered_top <- filtered %>%
  group_by(cluster, direction) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)), desc(abs(pct.1 - pct.2))) %>%
  slice_head(n = top) %>%
  ungroup()

# Save filtered markers
outfile <- file.path(outdir, paste0(id, "_markers_filtered.tsv"))
write.table(filtered_top, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)

# Summary table
write.table(
  data.frame(
    id,
    idents,
    n_cells_before,
    n_cells_after = ncol(obj),
    n_features_before,
    n_features_after = nrow(obj),
    dims,
    res,
    n_subclusters = length(unique(obj$seurat_clusters)),
    metric,
    seed,
    n_markers_total = nrow(markers),
    n_markers_filtered = nrow(filtered),
    n_markers_top = nrow(filtered_top),
    significance,
    regulation,
    enrichment,
    top
  ),
  file = file.path(outdir, paste0(id, "_subcluster_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
