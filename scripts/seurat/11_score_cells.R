#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(future)
})

set.seed(777)

# Set up parallelization
options(future.globals.maxSize = 16000 * 1024^2)  # 16 GB
plan("multicore", workers = 16)

# Close any open graphics devices
while (!is.null(dev.list())) dev.off()

# Args
args <- commandArgs(trailingOnly = TRUE)
id	<- args[1]
outdir	<- args[2]
input	<- args[3]

# Parameters
epi_thr <- 1.0
ery_thr <- 10

# Load clustered object
obj <- readRDS(input)

# Use SCT assay
DefaultAssay(obj) <- "SCT"

n_cells_before = ncol(obj)
n_features_before = nrow(obj)

# Define marker sets
epithelial_markers <- c("KRT8", "KRT18", "KRT1", "EPCAM", "CLDN7", 
                        "CDH1", "DSP", "DSG2", "EPPK1", "S100P")

erythrocyte_markers <- c("HBA1", "HBB", "HBA2", "HEMGN", 
                         "ALAS2", "PRDX2", "ANK1", "HBBP1")

schwann_markers <- c("SOX10", "S100", "S100B", "NGFR", "p75NTR", "MPZ", "MBP", "PMP22", "PLP1", "PRX", 
		     "NCAM", "NCAM1", "L1CAM", "SCN7A", "SOX2", "GAP43", "EGR2", "Krox20", "POU3F1", "OCT6")

epithelial_markers_available <- epithelial_markers[epithelial_markers %in% rownames(obj)]
erythrocyte_markers_available <- erythrocyte_markers[erythrocyte_markers %in% rownames(obj)]
schwann_markers_available <- schwann_markers[schwann_markers %in% rownames(obj)]

# Add module score for epithelial cells
obj <- AddModuleScore(
  obj,
  features = list(epithelial_markers_available),
  name = "Epithelial_score"
)
n_cells_epithelial = sum(obj$Epithelial_score1 >= epi_thr)

# Add module score for erythrocyte cells
obj <- AddModuleScore(
  obj,
  features = list(erythrocyte_markers_available),
  name = "Erythrocyte_score"
)
n_cells_erythrocyte = sum(obj$Erythrocyte_score1 >= ery_thr)

# Add module score for schwann cells
obj <- AddModuleScore(
  obj,
  features = list(schwann_markers_available),
  name = "Schwann_score"
)

# Plot VlnPlot
png(file.path(outdir, paste0(id, "_violin_score.png")), width = 1600, height = 1200)
VlnPlot(obj, features = c("Epithelial_score1", "Erythrocyte_score1", "Schwann_score1"), pt.size = 0, group.by = "orig.ident") +
  labs(x = NULL)
dev.off()

# Visualize distribution
png(file.path(outdir, paste0(id, "_epithelial_score_distribution.png")), width = 1600, height = 1200)
hist(obj$Epithelial_score1, breaks = 50, main = "Epithelial Score Distribution", xlab = "Epithelial Score")
abline(v = epi_thr, col = "red", lwd = 2, lty = 2)
dev.off()

png(file.path(outdir, paste0(id, "_erythrocyte_score_distribution.png")), width = 1600, height = 1200)
hist(obj$Erythrocyte_score1, breaks = 50, main = "Erythrocyte Score Distribution", xlab = "Erythrocyte Score")
abline(v = 10.0, col = "red", lwd = 2, lty = 2)
dev.off()

png(file.path(outdir, paste0(id, "_schwann_score_distribution.png")), width = 1600, height = 1200)
hist(obj$Schwann_score1, breaks = 50, main = "Schwann Score Distribution", xlab = "Schwann Score")
dev.off()

df <- obj@meta.data
png(file.path(outdir, paste0(id, "_scatter_plot.png")), width = 1600, height = 1200)
df %>% 
  arrange(Schwann_score1) %>%
  ggplot(aes_string("Epithelial_score1", "Erythrocyte_score1", colour = "Schwann_score1")) +
  geom_point() +
  scale_color_gradientn(colors = c("black","blue","green2","red","yellow"),
                        guide = guide_colorbar(barwidth = 2, barheight = 30)) +
  scale_x_log10() +
  scale_y_log10()
dev.off()

# Filtering
cells_to_keep <- which(obj$Epithelial_score1 < epi_thr & 
                       obj$Erythrocyte_score1 < ery_thr)

obj <- obj[, cells_to_keep]

# Summary table
write.table(
  data.frame(
    id,
    n_cells_before,
    n_cells_after = ncol(obj),
    n_features_before,
    n_features_after = nrow(obj),
    n_epithelial_markers_total = length(epithelial_markers),
    n_epithelial_markers_available = length(epithelial_markers_available),
    n_erythrocyte_markers_total = length(erythrocyte_markers),
    n_erythrocyte_markers_available = length(erythrocyte_markers_available),
    n_schwann_markers_total = length(schwann_markers),
    n_schwann_markers_available = length(schwann_markers_available),
    epi_thr,
    ery_thr,
    n_cells_epithelial,
    pct_cells_epithelial = n_cells_epithelial / n_cells_before * 100,
    n_cells_erythrocyte,
    n_cells_erythrocyte = n_cells_erythrocyte / n_cells_before * 100,
    pct_cells_kept = ncol(obj) / n_cells_before * 100
  ),
  file = file.path(outdir, paste0(id, "_enrichment_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)







# Parameters
npcs <- 50
dims <- 10
res <- 0.5
metric  <- "cosine"
seed    <- 777

significance <- 0.05
regulation <- 0.5
enrichment <- 0.2
top <- 30

# Use integrated assay
DefaultAssay(obj) <- "integrated"

# PCA
obj <- RunPCA(obj, npcs = npcs)

# Elbow plot
png(file.path(outdir, paste0(id, "_elbow.png")), width=1200, height=720)
print(ElbowPlot(obj, ndims = npcs) +
        geom_vline(xintercept = 10, linetype = "dashed", color = "red"))
dev.off()

# Clustering
obj <- FindNeighbors(obj, dims = 1:dims)
obj <- FindClusters(obj, res = res, algorithm = 4)

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
saveRDS(obj, file = file.path(outdir, paste0(id, "_enriched.rds")))

# Use SCT assay and prepare it
DefaultAssay(obj) <- "SCT"
obj <- PrepSCTFindMarkers(obj)

# Find markers
markers <- FindAllMarkers(
  obj,
  only.pos = FALSE
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
    n_cells = ncol(obj),
    n_features = nrow(obj),
    dims,
    res,
    n_clusters = length(unique(obj$seurat_clusters)),
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
  file = file.path(outdir, paste0(id, "_enriched_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
