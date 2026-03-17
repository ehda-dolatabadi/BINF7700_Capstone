#!/usr/bin/env Rscript
# Script: 14_auto_annotate.R
# Purpose: Annotate clusters
# Description: Annotate clusters automatically using SingleR
# Usage: Rscript 14_auto_annotate.R <id> <output> <input> <cells>

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
cells	<- unlist(strsplit(args[4], ","))

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
clusters <- obj@meta.data %>%
  group_by(seurat_clusters, singler_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(n, n = 1) %>%
  select(seurat_clusters, singler_label)

clusters <- setNames(clusters$singler_label, clusters$seurat_clusters)

# Add cluster identities as a metadata column
obj$singler_cluster <- unname(clusters[as.character(obj$seurat_clusters)])

# UMAP colored by SingleR cluster annotation
png(file.path(outdir, paste0(id, "_cluster_annotated.png")), width = 1600, height = 900, res=150)
print(DimPlot(obj, reduction = "umap", group.by = "singler_cluster", label = TRUE, repel = TRUE, label.size = 4) +
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

# Add cell counts to cell labels
cell_counts <- table(obj$singler_label)
cell_labels <- paste0(names(cell_counts), " (", cell_counts, ")")
names(cell_labels) <- names(cell_counts)

# UMAP colored by SingleR cell annotation
png(file.path(outdir, paste0(id, "_cell_annotated.png")), width = 1600, height = 900, res=150)
print(DimPlot(obj, reduction = "umap", group.by = "singler_label", label = FALSE) +
	scale_color_discrete(labels = cell_labels) +
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

# Filter desired cell type
cells_to_keep <- split(rownames(obj@meta.data), obj$singler_label) [cells]

# Filtered UMAP colored by SingleR cell annotation
png(file.path(outdir, paste0(id, "_cell_annotated_filtered.png")), width = 1600, height = 900, res=150)
print(DimPlot(obj, reduction = "umap", cells.highlight = cells_to_keep, cols.highlight = scales::hue_pal()(length(cells)), cols = "lightgrey") + 
        labs(
                title = "UMAP Clustering",
                x = "UMAP 1",
                y = "UMAP 2"
        ) +
	guides(color = guide_legend(ncol = 1)) +
        theme(
                plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                axis.title = element_text(size = 14, face = "bold"),
                axis.text = element_text(size = 12),
		legend.key.width = unit(1, "cm"),
                legend.text = element_text(size = 11)
        ))
dev.off()

# Save object
saveRDS(obj, file = file.path(dirname(outdir), paste0(id, "_processed.rds")))
