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

# Join split sample layers into one counts matrix (each sample is stored separately)
obj <- JoinLayers(obj, assay = "RNA")

# Set default assay to RNA to ensure full gene set is transferred
DefaultAssay(obj) <- "RNA"

# Convert Seurat object to Monocle3 cell_data_set
# this transfers raw counts, cell metadata, and gene metadata
cds <- as.cell_data_set(obj)

# Normalize and reduce to 10 dimensions using PCA
# num_dim = 10 based on meaningful components identified in prior Seurat analysis
cds <- preprocess_cds(cds, num_dim = 10)

# Compute UMAP embedding using Monocle3's own pipeline
# this is separate from Seurat's UMAP and optimized for trajectory inference
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")

# Cluster cells using Louvain community detection on the UMAP embedding
# Monocle3 uses these clusters internally to structure the trajectory graph
cds <- cluster_cells(cds, reduction_method = "UMAP")

# Learn the principal graph separately for each partition
# use_partition = TRUE allows disconnected trajectories for unrelated cell populations
cds <- learn_graph(cds, use_partition = FALSE)

# Select control cells as the pseudotime root based on orig.ident metadata
earliest_cells <- rownames(colData(cds)[colData(cds)$orig.ident == "control", ])

# Assign pseudotime values to all cells starting from the control cells
cds <- order_cells(cds, root_cells = earliest_cells)

# Plot trajectory colored by pseudotime value
png(file.path(outdir, paste0(id, "_trajectory_pseudotime.png")), width = 1600, height = 1200, res = 150)
print(plot_cells(cds,
                 color_cells_by      = "pseudotime",     # color by pseudotime value
                 label_cell_groups   = FALSE,            # suppress cluster ID labels on plot
                 label_leaves        = TRUE,             # label trajectory endpoints
                 label_branch_points = FALSE))            # label trajectory branch points
dev.off()

# Save full CDS object including nearest neighbor indices
# save_monocle_objects is preferred over saveRDS to preserve all internal graph structures
save_monocle_objects(cds            = cds,
                     directory_path = file.path(outdir, paste0(id, "_trajectory")))
