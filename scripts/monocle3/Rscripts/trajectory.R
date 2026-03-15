#!/usr/bin/env Rscript
# Script: trajectory.R
# Purpose:
# Description:
# Usage: Rscript trajectory.R <id> <output> <input>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(remotes)
})

if (!requireNamespace("SeuratWrappers", quietly = TRUE))
    remotes::install_github("satijalab/seurat-wrappers")

suppressPackageStartupMessages({
  library(SeuratWrappers)
})

set.seed(777)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]

# Load object
obj <- readRDS(input)

# Convert to CDS
cds <- as.cell_data_set(obj)

# Recreate partitions (all in one partition for connected trajectory)
partitions <- rep(1, ncol(cds))
names(partitions) <- colnames(cds)
cds@clusters$UMAP$partitions <- as.factor(partitions)

# Transfer cluster labels from Seurat
cds@clusters$UMAP$clusters <- Idents(obj)

# Transfer UMAP coordinates
cds@int_colData@listData$reducedDims$UMAP <- obj@reductions$umap@cell.embeddings

# Learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

# Order cells by pseudotime
earliest_cells <- rownames(colData(cds)[colData(cds)$orig.ident == "control", ])
cds <- order_cells(cds, root_cells = earliest_cells)

# Plot by cell type
png(file.path(outdir, paste0(id, "_trajectory_celltype.png")), width = 1600, height = 1200, res=150)
print(plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster = TRUE,
           label_leaves = TRUE,
           label_branch_points = TRUE))
dev.off()

# Plot by pseudotime
png(file.path(outdir, paste0(id, "_trajectory_pseudotime.png")), width = 1600, height = 1200, res=150)
print(plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE))
dev.off()

# Save object
saveRDS(cds, file = file.path(outdir, paste0(id, "_trajectory.rds")))
