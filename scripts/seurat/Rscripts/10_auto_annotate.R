#!/usr/bin/env Rscript
# Script: 10_auto_annotate.R
# Purpose: Annotate clusters
# Description: Annotate clusters automatically using SingleR
# Usage: Rscript 10_auto_annotate.R <id> <output> <input> <cells>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(SingleR)
  library(celldex)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(271)

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
# SingleR has no seed parameter; uses bootstrap internally — set.seed() here influences bootstrap sampling
set.seed(271)
predictions <- SingleR(
  test      = LayerData(obj, assay = "SCT", layer = "data"),
  ref       = ref,
  labels    = ref$label.main,
  fine.tune = TRUE
)

# Add labels to Seurat object
obj$singler_label	<- predictions$labels
obj$singler_pruned	<- predictions$pruned.labels
obj$singler_score	<- apply(predictions$scores, 1, max)

# plots
for (label_type in c("singler_label", "singler_pruned")) {

  suffix      <- ifelse(label_type == "singler_label", "label", "pruned")
  cluster_col <- paste0("singler_cluster_", suffix)

  # Get consensus label per cluster by majority vote
  # NA cells (pruned) are excluded from the vote
  clusters <- obj@meta.data %>%
    filter(!is.na(.data[[label_type]])) %>%
    group_by(seurat_clusters, .data[[label_type]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    select(seurat_clusters, all_of(label_type))
  clusters <- setNames(clusters[[label_type]], clusters$seurat_clusters)

  # Add cluster identities as a metadata column
  obj[[cluster_col]] <- unname(clusters[as.character(obj$seurat_clusters)])

  # UMAP colored by SingleR cluster annotation
  png(file.path(outdir, paste0(id, "_cluster_annotated_", suffix, ".png")), width = 1600, height = 900, res = 150)
  print(DimPlot(obj, reduction = "umap", group.by = cluster_col, label = TRUE, repel = TRUE, label.size = 4) +
    labs(
      title = "UMAP Clustering",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme(
      plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title       = element_text(size = 14, face = "bold"),
      axis.text        = element_text(size = 12),
      legend.key.width = unit(1.5, "cm"),
      legend.text      = element_text(size = 11)
    ))
  dev.off()

  # Add cell counts to cell labels (table excludes NA by default)
  cell_counts <- table(obj@meta.data[[label_type]])
  cell_labels <- paste0(names(cell_counts), " (", cell_counts, ")")
  names(cell_labels) <- names(cell_counts)

  # Order labels by cell count (descending)
  label_counts <- sort(table(obj@meta.data[[label_type]]), decreasing = TRUE)
  obj@meta.data[[label_type]] <- factor(
    obj@meta.data[[label_type]],
    levels = names(label_counts)
  )

  # Rebuild cell_labels in the same order
  cell_labels <- paste0(names(label_counts), " (", label_counts, ")")
  names(cell_labels) <- names(label_counts)

  # UMAP colored by SingleR cell annotation
  png(file.path(outdir, paste0(id, "_cell_annotated_", suffix, ".png")), width = 1600, height = 900, res = 150)
  print(DimPlot(obj, reduction = "umap", group.by = label_type, label = FALSE) +
    scale_color_discrete(labels = cell_labels, na.value = "lightgrey") +
    labs(
      title = "UMAP Clustering",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme(
      plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title       = element_text(size = 14, face = "bold"),
      axis.text        = element_text(size = 12),
      legend.key.width = unit(1.5, "cm"),
      legend.text      = element_text(size = 8)
    ))
  dev.off()
}

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
    theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
		legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size = 11)
    ))
dev.off()

# Save object
saveRDS(obj, file = file.path(dirname(outdir), paste0(id, "_processed.rds")))
