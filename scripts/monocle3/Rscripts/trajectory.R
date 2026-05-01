#!/usr/bin/env Rscript
# Script: trajectory.R
# Purpose: Infer pseudotime trajectories from Monocle3 for axolotl limb regeneration scRNA-seq data
# Description: Runs a separate full Monocle3 pipeline per lineage (defined by singler_cluster),
#              so each lineage gets its own UMAP, graph, and pseudotime scale. Saves diagnostic
#              plots, graph test results, and model per lineage. Produces one summary row per lineage.
# Usage: Rscript trajectory.R <id> <output> <input>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(ggplot2)
  library(remotes)
})

if (!requireNamespace("SeuratWrappers", quietly = TRUE))
    remotes::install_github("satijalab/seurat-wrappers")

suppressPackageStartupMessages({
  library(SeuratWrappers)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(271)

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

# Define lineages as named groups of singler_cluster labels
# each lineage gets an independent Monocle3 run so its UMAP reflects its own structure
# rather than being distorted by the dominant myeloid population
lineages <- list(
  myeloid    = c("Monocyte", "DC", "Pro-Myelocyte"),
  neutrophil = c("Neutrophils"),
  tcell      = c("T_cells"),
  neuronal   = c("Neurons")
)

# Ordered timepoints used for per-lineage root selection
timepoint_order <- c("control", "3h", "24h", "72h", "7dpa", "14dpa", "22dpa", "33dpa")

# Accumulate one summary row per lineage
summary_rows <- list()

# =============================================================================
# Per-lineage Monocle3 pipeline
# =============================================================================

for (lineage_name in names(lineages)) {

  cluster_labels <- lineages[[lineage_name]]
  lin_id         <- paste0(id, "_", lineage_name)
  lin_outdir     <- file.path(outdir, lineage_name)
  dir.create(lin_outdir, showWarnings = FALSE, recursive = TRUE)

  # Subset Seurat object to cells belonging to this lineage
  obj_lin <- subset(obj, subset = singler_cluster %in% cluster_labels)

  # Skip lineage if too few cells for meaningful trajectory inference
  if (ncol(obj_lin) < 50) next

  # Convert lineage subset to Monocle3 cell_data_set
  # each lineage is converted independently so preprocessing reflects its own HVGs
  cds <- as.cell_data_set(obj_lin)

  # Normalize and reduce to 50 dimensions using PCA
  # preprocess_cds has no seed parameter; calls set.seed(2016) internally — external set.seed() is overridden and ineffective
  set.seed(271)
  cds <- preprocess_cds(cds, num_dim = 50)

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

  # ===========================================================================
  # Diagnostic plots
  # ===========================================================================

  # Scale dot size inversely with cell count so small lineages remain visible;
  # capped between 0.5 (large populations) and the natural scaling for small ones
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

  # Helper to build per-group scale labels with cell counts
  make_labels <- function(col) {
    counts <- table(col)
    labels <- paste0(names(counts), " (", counts, ")")
    setNames(labels, names(counts))
  }

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
    scale_color_discrete(labels = make_labels(colData(cds)$orig.ident),
                         na.value = "lightgrey") +
    theme(aspect.ratio = 0.75)))
  dev.off()

  # Plot 3: Cluster-level singler annotation
  png(file.path(lin_outdir, paste0(lin_id, "_trajectory_singler_cluster.png")),
      width = 1600, height = 900, res = 150)
  print(remove_black_outline(
    do.call(plot_cells, c(list(cds = cds, color_cells_by = "singler_cluster"), plot_args)) +
    scale_color_discrete(labels = make_labels(colData(cds)$singler_cluster),
                         na.value = "lightgrey") +
    theme(aspect.ratio = 0.75)))
  dev.off()

  # Plot 4: Cell-level singler annotation
  png(file.path(lin_outdir, paste0(lin_id, "_trajectory_singler_label.png")),
      width = 1600, height = 900, res = 150)
  print(remove_black_outline(
    do.call(plot_cells, c(list(cds = cds, color_cells_by = "singler_label"), plot_args)) +
    scale_color_discrete(labels = make_labels(colData(cds)$singler_label),
                         na.value = "lightgrey") +
    theme(aspect.ratio = 0.75)))
  dev.off()

  # Plot 5: Seurat cluster numbers with numeric sort to prevent lexicographic ordering
  colData(cds)$seurat_clusters <- factor(
    colData(cds)$seurat_clusters,
    levels = as.character(sort(as.integer(levels(colData(cds)$seurat_clusters))))
  )
  png(file.path(lin_outdir, paste0(lin_id, "_trajectory_seurat_clusters.png")),
      width = 1600, height = 900, res = 150)
  print(remove_black_outline(
    do.call(plot_cells, c(list(cds = cds, color_cells_by = "seurat_clusters"), plot_args)) +
    scale_color_discrete(labels = make_labels(colData(cds)$seurat_clusters),
                         na.value = "lightgrey") +
    theme(aspect.ratio = 0.75)))
  dev.off()

  # ===========================================================================
  # Trajectory-variable gene analysis
  # ===========================================================================

  # Test for genes whose expression varies significantly along the principal graph
  # using Moran's I spatial autocorrelation on the graph topology
  # graph_test is stochastic when cores > 1 due to non-deterministic parallel scheduling;
  # results may vary slightly across runs at the gene ranking margins
  set.seed(271)
  graph_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

  write.table(
    graph_test_res[order(graph_test_res$q_value), ],
    file      = file.path(lin_outdir, paste0(lin_id, "_trajectory_graph_test.tsv")),
    sep       = "\t",
    quote     = FALSE,
    row.names = TRUE
  )

  traj_deg_ids <- rownames(subset(graph_test_res, q_value < 0.05 & morans_I > 0.1))

  # ===========================================================================
  # Save model
  # ===========================================================================

  # save_monocle_objects is preferred over saveRDS to preserve all internal graph structures
  save_monocle_objects(cds            = cds,
                       directory_path = file.path(lin_outdir,
                                                  paste0(lin_id, "_trajectory")))

  # ===========================================================================
  # Accumulate summary row for this lineage
  # ===========================================================================

  timepoint_counts       <- table(colData(cds)$orig.ident)
  singler_cluster_counts <- table(colData(cds)$singler_cluster)

  summary_rows[[lineage_name]] <- data.frame(
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
  )
}

# =============================================================================
# Summary
# =============================================================================

write.table(
  do.call(rbind, summary_rows),
  file      = file.path(outdir, paste0(id, "_trajectory_summary.tsv")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)