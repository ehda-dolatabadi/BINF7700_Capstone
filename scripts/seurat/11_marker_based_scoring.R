#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]
dim	<- args[4]
res	<- args[5]

# Load clustered object
obj <- readRDS(input)

# Core Schwann cell identity markers
schwann_markers <- c("PRX", "MBP", "MPZ", "SOX10", "PLP1", "S100B")
schwann_markers_available <- schwann_markers[schwann_markers %in% rownames(obj)]

# Add module score (creates "Schwann1" in metadata)
obj <- AddModuleScore(
  obj,
  features = list(schwann_markers_available),
  name = "Schwann"
)

# Visualize scores
FeaturePlot(obj, features = "Schwann1")
VlnPlot(obj, features = "Schwann1", group.by = "seurat_clusters")

















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

# Re-run clustering
obj <- FindNeighbors(obj, dims = 1:dims)
obj <- FindClusters(obj, res = res, algorithm = 4)

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_clustered.rds")))

# Summary
write.table(
  data.frame(
    id		= id,
    n_cells	= ncol(obj),
    n_features	= nrow(obj),
    dims_used	= dims,
    resolution	= res,
    n_clusters	= length(levels(obj$seurat_clusters))
  ),
  file = file.path(outdir, paste0(id, "_clustered.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
