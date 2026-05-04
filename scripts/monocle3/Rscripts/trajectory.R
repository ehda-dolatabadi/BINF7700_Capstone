#!/usr/bin/env Rscript
# Script: trajectory.R
# Purpose: Infer pseudotime trajectories from Monocle3 for axolotl limb regeneration scRNA-seq data
# Description: Runs a full Monocle3 pipeline for one lineage defined by singler_cluster labels,
#              so each lineage gets its own UMAP, graph, and pseudotime scale. Saves diagnostic
#              plots, graph test results, and model. Called once per lineage via SLURM array job.
# Usage: Rscript trajectory.R <id> <output> <input> <lineage> <clusters> <dims>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(ggplot2)
  library(SeuratWrappers)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(271)

# Parse command line arguments
args           <- commandArgs(trailingOnly = TRUE)
id             <- args[1]
outdir         <- args[2]
input          <- args[3]
lineage_name   <- args[4]                        # e.g. "myeloid"
cluster_labels <- strsplit(args[5], ";")[[1]]    # e.g. "Monocyte;DC;Pro-Myelocyte"
num_dim        <- as.numeric(args[6])

# Ordered timepoints used for root selection
timepoint_order <- c("control", "3h", "24h", "72h", "7dpa", "14dpa", "22dpa", "33dpa")

# Per-lineage output identifiers and directories
lin_id     <- paste0(id, "_", lineage_name)
lin_outdir <- file.path(outdir, lineage_name)
dir.create(lin_outdir, showWarnings = FALSE, recursive = TRUE)

# Load object
obj <- readRDS(input)

# Join split sample layers into one counts matrix (each sample is stored separately)
obj <- JoinLayers(obj, assay = "RNA")

# Set default assay to RNA to ensure full gene set is transferred
DefaultAssay(obj) <- "RNA"

# Subset Seurat object to cells belonging to this lineage
obj_lin <- subset(obj, subset = singler_cluster %in% cluster_labels)

# Convert lineage subset to Monocle3 cell_data_set
# each lineage is converted independently so preprocessing reflects its own HVGs
cds <- as.cell_data_set(obj_lin)

# Normalize and reduce to num_dim dimensions using PCA
# preprocess_cds has no seed parameter; calls set.seed(2016) internally — external set.seed() is overridden and ineffective
set.seed(271)
cds <- preprocess_cds(cds, num_dim = num_dim)

# Compute UMAP embedding for this lineage
# reduce_dimension has no seed parameter; calls set.seed(2016) internally — external set.seed() is overridden and ineffective
set.seed(271)
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")

# Cluster cells using Louvain community detection on the lineage UMAP
cds <- cluster_cells(cds, reduction_method = "UMAP", random_seed = 271)  # default: 42

# Learn principal graph within this lineage
# use_partition = FALSE: biological separation was already done by subsetting;
# forcing a single connected graph within each lineage is appropriate here
# learn_graph has no seed parameter; uses internal graph learning with stochastic initialization
set.seed(271)
cds <- learn_graph(cds, use_partition = FALSE)

# Select root cells from the earliest timepoint present in this lineage subset
cell_meta <- as.data.frame(colData(cds))
earliest_cells <- character(0)
for (tp in timepoint_order) {
  candidates <- rownames(cell_meta[cell_meta$orig.ident == tp, ])
  if (length(candidates) > 0) {
    earliest_cells <- candidates
    break
  }
}

# Assign pseudotime from the earliest timepoint root
cds <- order_cells(cds, root_cells = earliest_cells)

# Reorder a colData column as a factor sorted by count descending
# so ggplot legend order reflects cell abundance
order_by_count <- function(col) {
  counts <- sort(table(col), decreasing = TRUE)
  factor(col, levels = names(counts))
}

# Helper to build per-group scale labels with cell counts
make_labels <- function(col) {
  counts <- sort(table(col), decreasing = TRUE)
  labels <- paste0(names(counts), " (", counts, ")")
  setNames(labels, names(counts))
}

# Remove the hardcoded black background geom_point layer from plot_cells output;
# ggplot2 stores color internally as "colour" not "color"
remove_black_outline <- function(p) {
  p$layers <- p$layers[!sapply(p$layers, function(l) {
    inherits(l$geom, "GeomPoint") &&
    !is.null(l$aes_params$colour) &&
    l$aes_params$colour == "black"
  })]
  p
}

# Generate a gradient color palette for n levels to avoid similar hues from
# the default discrete palette when there are many categories
gradient_colors <- function(levels) {
  setNames(
    colorRampPalette(c("darkblue", "purple", "magenta", "orange", "yellow"))(length(levels)),
    levels
  )
}

# Reorder colData columns for consistent legend ordering
colData(cds)$orig.ident      <- factor(colData(cds)$orig.ident, levels = timepoint_order)
colData(cds)$singler_cluster <- order_by_count(factor(colData(cds)$singler_cluster))
colData(cds)$singler_label   <- order_by_count(factor(colData(cds)$singler_label))
colData(cds)$seurat_clusters <- factor(
  colData(cds)$seurat_clusters,
  levels = as.character(sort(as.integer(levels(colData(cds)$seurat_clusters))))
)

# Scale dot size linearly from 1 (at 6000+ cells) to 3 (at 0 cells)
cell_size <- max(1, min(3, 3 - (ncol(cds) / 3000)))

# Shared plot arguments across all plots
# cell_stroke = I(0) removes point border; the black outline is from a separate
# black geom_point background layer hardcoded in plot_cells at 1.5x cell_size —
# it cannot be disabled via parameters and must be removed from the ggplot layers
plot_args <- list(
  label_cell_groups       = FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves            = FALSE,
  label_branch_points     = FALSE,
  label_roots             = FALSE,
  label_principal_points  = FALSE,
  show_trajectory_graph   = TRUE,
  cell_size               = cell_size,
  cell_stroke             = I(0)
)

# Plot 1: Pseudotime with trajectory graph — reference plot for matching node
# positions against all other plots; graph lines shown on all plots
png(file.path(lin_outdir, paste0(lin_id, "_trajectory_pseudotime.png")),
    width = 1600, height = 900, res = 150)
print(remove_black_outline(
  do.call(plot_cells, c(list(cds = cds, color_cells_by = "pseudotime"), plot_args)) +
  theme(aspect.ratio = 0.75)))
dev.off()

# Plot 2: Timepoint
png(file.path(lin_outdir, paste0(lin_id, "_trajectory_timepoint.png")),
    width = 1600, height = 900, res = 150)
print(remove_black_outline(
  do.call(plot_cells, c(list(cds = cds, color_cells_by = "orig.ident"), plot_args)) +
  scale_color_manual(
    values = gradient_colors(levels(colData(cds)$orig.ident)),
    labels = make_labels(colData(cds)$orig.ident),
    na.value = "lightgrey"
  ) +
  theme(aspect.ratio = 0.75)))
dev.off()

# Plot 3: Cluster-level singler annotation
png(file.path(lin_outdir, paste0(lin_id, "_trajectory_singler_cluster.png")),
    width = 1600, height = 900, res = 150)
print(remove_black_outline(
  do.call(plot_cells, c(list(cds = cds, color_cells_by = "singler_cluster"), plot_args)) +
  scale_color_manual(
    values = gradient_colors(levels(colData(cds)$singler_cluster)),
    labels = make_labels(colData(cds)$singler_cluster),
    na.value = "lightgrey"
  ) +
  theme(aspect.ratio = 0.75)))
dev.off()

# Plot 4: Cell-level singler annotation
png(file.path(lin_outdir, paste0(lin_id, "_trajectory_singler_label.png")),
    width = 1600, height = 900, res = 150)
print(remove_black_outline(
  do.call(plot_cells, c(list(cds = cds, color_cells_by = "singler_label"), plot_args)) +
  scale_color_manual(
    values = gradient_colors(levels(colData(cds)$singler_label)),
    labels = make_labels(colData(cds)$singler_label),
    na.value = "lightgrey"
  ) +
  theme(aspect.ratio = 0.75)))
dev.off()

# Plot 5: Seurat cluster numbers with numeric sort to prevent lexicographic ordering
png(file.path(lin_outdir, paste0(lin_id, "_trajectory_seurat_clusters.png")),
    width = 1600, height = 900, res = 150)
print(remove_black_outline(
  do.call(plot_cells, c(list(cds = cds, color_cells_by = "seurat_clusters"), plot_args)) +
  scale_color_manual(
    values = gradient_colors(levels(colData(cds)$seurat_clusters)),
    labels = make_labels(colData(cds)$seurat_clusters),
    na.value = "lightgrey"
  ) +
  theme(aspect.ratio = 0.75)))
dev.off()

# Test for genes whose expression varies significantly along the principal graph
# using Moran's I spatial autocorrelation on the graph topology
# graph_test is stochastic when cores > 1 due to non-deterministic parallel scheduling;
# results may vary slightly across runs at the gene ranking margins
set.seed(271)
graph_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)

write.table(
  graph_test_res[order(graph_test_res$q_value), ],
  file      = file.path(lin_outdir, paste0(lin_id, "_trajectory_graph_test.tsv")),
  sep       = "\t",
  quote     = FALSE,
  row.names = TRUE
)

traj_deg_ids <- rownames(subset(graph_test_res, q_value < 0.05 & morans_I > 0.1))

# save_monocle_objects is preferred over saveRDS to preserve all internal graph structures
save_monocle_objects(cds            = cds,
                     directory_path = file.path(lin_outdir,
                                                paste0(lin_id, "_trajectory")))

timepoint_counts       <- table(colData(cds)$orig.ident)
singler_cluster_counts <- table(colData(cds)$singler_cluster)

write.table(
  data.frame(
    id                        = lin_id,
    lineage                   = lineage_name,
    n_cells                   = ncol(cds),
    n_monocle_clusters        = length(unique(clusters(cds))),
    n_partitions              = length(unique(partitions(cds))),
    root_timepoint            = timepoint_order[min(which(timepoint_order %in%
                                  colData(cds)$orig.ident))],
    n_root_cells              = length(earliest_cells),
    n_trajectory_genes        = length(traj_deg_ids),
    cells_per_timepoint       = paste(timepoint_order,
                                      as.integer(timepoint_counts[timepoint_order]),
                                      sep = ":", collapse = "; "),
    cells_per_singler_cluster = paste(names(singler_cluster_counts),
                                      singler_cluster_counts,
                                      sep = ":", collapse = "; ")
  ),
  file      = file.path(lin_outdir, paste0(lin_id, "_trajectory_summary.tsv")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)