#!/usr/bin/env Rscript
# Script: 14_auto_annotate.R
# Purpose: Annotate clusters
# Description: Annotate clusters automatically using SingleR
# Usage: Rscript 14_auto_annotate.R <id> <output> <input>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(SingleR)
  library(celldex)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]

# Load object
obj <- readRDS(input)

# Load reference dataset
#	HumanPrimaryCellAtlasData()		broad human cell types
#	BlueprintEncodeData()			human, good for immune and stromal
#	MouseRNAseqData()			mouse bulk RNA reference
#	ImmGenData()				mouse immune cells (high resolution)
#	DatabaseImmuneCellExpressionData()	human immune cells

ref <- celldex::HumanPrimaryCellAtlasData()

# Run SingleR
#	label.main: broad cell types
#	label.fine: finer subtypes
#	fine.tune: TRUE improves accuracy on ambiguous cells

predictions <- SingleR(
  test      = obj@assays$SCT@data,
  ref       = ref,
  labels    = ref$label.main,
  fine.tune = TRUE
)

# Add labels to Seurat object
obj$singler_label	<- predictions$labels
obj$singler_pruned	<- predictions$pruned.labels
obj$singler_score	<- apply(predictions$scores, 1, max)

# Get consensus label per cluster by majority vote
cluster_labels <- obj@meta.data %>%
  group_by(seurat_clusters, singler_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(n, n = 1) %>%
  select(seurat_clusters, singler_label)

cluster_labels <- setNames(cluster_labels$singler_label, cluster_labels$seurat_clusters)

# Rename cluster identities in Seurat object
obj <- RenameIdents(obj, cluster_labels)

# Plots
png(file.path(outdir, paste0(id, "_auto_annotated.png")), width = 1600, height = 1200, res=150)
print(DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5) +
        labs(
                title = "UMAP Clustering",
                x = "UMAP 1",
                y = "UMAP 2"
        ) +
        theme(
                plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                axis.title = element_text(size = 14, face = "bold"),
                axis.text = element_text(size = 12),
                legend.text = element_text(size = 11)
        ))
dev.off()

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_auto_annotated.rds")))
