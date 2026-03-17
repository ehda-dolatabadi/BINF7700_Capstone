#!/usr/bin/env Rscript
# Script: 13_manual_annotate.R
# Purpose: Annotate clusters
# Description: Renames cluster identities using ordered labels
# Usage: Rscript 13_manual_annotate.R <id> <output> <input> <labels>
#	<labels>: comma-separated labels in cluster order

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]
labels	<- unlist(strsplit(args[4], ","))
labels	<- setNames(labels, seq_along(labels))

# Load object
obj <- readRDS(input)

# Annotate
obj <- RenameIdents(obj, labels)

# Plot
png(file.path(outdir, paste0(id, "_annotated.png")), width = 1600, height = 1200, res=150)
print(DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4) +
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
saveRDS(obj, file = file.path(dirname(outdir), paste0(id, "_processed.rds")))
