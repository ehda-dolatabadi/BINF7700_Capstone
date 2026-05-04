#!/usr/bin/env Rscript
# Script: cellchat.R
# Purpose: Infer ligand-receptor mediated cell-cell communication networks per timepoint
# Description: Loads an annotated Seurat object, splits it by timepoint (orig.ident), and
#              for each timepoint runs the standard CellChat v2 inference pipeline
#              (over-expressed genes/interactions, communication probability, pathway
#              aggregation, network centrality). Outputs per-timepoint CellChat objects
#              and plots, plus a merged object with cross-timepoint comparison plots.
# Usage: Rscript cellchat.R <id> <outdir> <input> <group_by>
#   <group_by>: metadata column to use as cell groups (e.g. "singler_cluster")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(remotes)
})
if (!requireNamespace("NMF", quietly = TRUE))
    remotes::install_version("NMF", version = "0.28.0", repos = "https://cloud.r-project.org", upgrade = "never")
if (!requireNamespace("CellChat", quietly = TRUE))
    remotes::install_github("jinworks/CellChat@v2.1.2", upgrade = "never")

suppressPackageStartupMessages({
  library(CellChat)
})
if (!requireNamespace("presto", quietly = TRUE)) {
    Sys.setenv("PKG_CXXFLAGS" = "-std=gnu++14")
    remotes::install_github("immunogenomics/presto@1.0.0", upgrade = "never")
}

RNGkind("L'Ecuyer-CMRG")
set.seed(271)

# Configure parallel processing
options(future.seed = TRUE)
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 16)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id       <- args[1]
outdir   <- args[2]
input    <- args[3]
group_by <- args[4]

# Load annotated Seurat object
obj <- readRDS(input)
DefaultAssay(obj) <- "SCT"

# Timepoint order from aggr_samples.csv
timepoint_order <- c("control", "3h", "24h", "72h", "7dpa", "14dpa", "22dpa", "33dpa")
timepoints <- intersect(timepoint_order, unique(as.character(obj$orig.ident)))

# Per-timepoint CellChat inference
cellchat_list <- list()
summary_rows  <- list()

for (tp in timepoints) {

  tp_id     <- paste0(id, "_", tp)
  tp_outdir <- file.path(outdir, tp_id)
  dir.create(tp_outdir, recursive = TRUE, showWarnings = FALSE)

  # Subset to this timepoint
  obj_tp <- subset(obj, subset = orig.ident == tp)
  if (ncol(obj_tp) < 100) {
    message("Skipping ", tp, ": only ", ncol(obj_tp), " cells")
    next
  }

  # Create CellChat object using human ligand-receptor database
  # (axolotl symbols are mapped to human via loc_map.tsv)
  cellchat <- createCellChat(object = obj_tp, group.by = group_by, assay = "SCT")
  cellchat@DB <- CellChatDB.human

  # Identify over-expressed signaling genes and interactions
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # Infer communication probabilities at L-R and pathway level
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)

  # Aggregate communication network and compute centrality scores
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  # Aggregated interaction count circle plot
  png(file.path(tp_outdir, paste0(tp_id, "_interactions_count.png")),
      width = 1200, height = 1200, res = 150)
  netVisual_circle(
    cellchat@net$count,
    vertex.weight = as.numeric(table(cellchat@idents)),
    weight.scale  = TRUE,
    label.edge    = FALSE,
    title.name    = paste0(tp, ": Number of interactions")
  )
  dev.off()

  # Aggregated interaction strength circle plot
  png(file.path(tp_outdir, paste0(tp_id, "_interactions_weight.png")),
      width = 1200, height = 1200, res = 150)
  netVisual_circle(
    cellchat@net$weight,
    vertex.weight = as.numeric(table(cellchat@idents)),
    weight.scale  = TRUE,
    label.edge    = FALSE,
    title.name    = paste0(tp, ": Interaction strength")
  )
  dev.off()

  # Signaling role scatter: outgoing vs incoming per cell group
  png(file.path(tp_outdir, paste0(tp_id, "_signaling_role.png")),
      width = 1600, height = 1200, res = 150)
  print(netAnalysis_signalingRole_scatter(cellchat) +
        ggtitle(tp) +
        theme(aspect.ratio = 0.75))
  dev.off()

  # Pathway-by-cell-group heatmaps for outgoing and incoming signaling
  png(file.path(tp_outdir, paste0(tp_id, "_pathway_heatmap_outgoing.png")),
      width = 1200, height = 1800, res = 150)
  print(netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing"))
  dev.off()

  png(file.path(tp_outdir, paste0(tp_id, "_pathway_heatmap_incoming.png")),
      width = 1200, height = 1800, res = 150)
  print(netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming"))
  dev.off()

  # Save per-timepoint CellChat object
  saveRDS(cellchat, file = file.path(tp_outdir, paste0(tp_id, "_cellchat.rds")))
  cellchat_list[[tp]] <- cellchat

  # Accumulate summary row
  summary_rows[[tp]] <- data.frame(
    id                   = tp_id,
    timepoint            = tp,
    group_by             = group_by,
    n_cells              = length(cellchat@idents),
    n_groups             = length(levels(cellchat@idents)),
    n_signaling_genes    = length(cellchat@var.features$features),
    n_LR_pairs           = nrow(cellchat@LR$LRsig),
    n_significant_LR     = sum(cellchat@net$pval < 0.05, na.rm = TRUE),
    n_signaling_pathways = length(cellchat@netP$pathways)
  )
}

# Cross-timepoint comparison
if (length(cellchat_list) >= 2) {

  cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

  # Total interaction count across timepoints
  png(file.path(outdir, paste0(id, "_compare_interactions_count.png")),
      width = 1600, height = 1200, res = 150)
  print(compareInteractions(cellchat_merged, show.legend = FALSE,
                            group = seq_along(cellchat_list), measure = "count"))
  dev.off()

  # Total interaction strength across timepoints
  png(file.path(outdir, paste0(id, "_compare_interactions_weight.png")),
      width = 1600, height = 1200, res = 150)
  print(compareInteractions(cellchat_merged, show.legend = FALSE,
                            group = seq_along(cellchat_list), measure = "weight"))
  dev.off()

  # Pathway information flow ranking across timepoints
  png(file.path(outdir, paste0(id, "_rank_net.png")),
      width = 1200, height = 1800, res = 150)
  print(rankNet(cellchat_merged, mode = "comparison",
                stacked = TRUE, do.stat = TRUE))
  dev.off()

  # Save merged object
  saveRDS(cellchat_merged,
          file = file.path(outdir, paste0(id, "_cellchat_merged.rds")))
}

# Summary table
write.table(
  do.call(rbind, summary_rows),
  file      = file.path(outdir, paste0(id, "_cellchat_summary.tsv")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
