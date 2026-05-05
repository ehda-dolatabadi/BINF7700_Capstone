#!/usr/bin/env Rscript
# Script: inference_per_timepoint.R
# Purpose: Infer ligand-receptor mediated cell-cell communication for one timepoint
# Description: Loads an annotated Seurat object, subsets to the specified timepoint,
#              and runs the standard CellChat v2 inference pipeline (over-expressed
#              genes/interactions, communication probability, pathway aggregation,
#              network centrality). Saves the CellChat object, plots, and a one-row
#              summary TSV. Designed to run as one task in an sbatch array job.
# Usage: Rscript inference_per_timepoint.R <id> <outdir> <input> <group_by> <tp>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(future)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(271)

# Configure parallel processing
options(future.seed = TRUE)
options(future.globals.maxSize = 48000 * 1024^2)  # 48 GB
plan("multicore", workers = 16)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id       <- args[1]
outdir   <- args[2]
input    <- args[3]
group_by <- args[4]
tp       <- args[5]

tp_id <- paste0(id, "_", tp)

# Load annotated Seurat object and subset to this timepoint
obj <- readRDS(input)
DefaultAssay(obj) <- "SCT"

obj <- subset(obj, subset = orig.ident == tp)
if (ncol(obj) < 100) {
  message("Skipping ", tp, ": only ", ncol(obj), " cells")
  quit(save = "no", status = 0)
}

# Drop cell-group levels not present in this timepoint
# (CellChat fails on empty factor levels)
obj@meta.data[[group_by]] <- droplevels(factor(obj@meta.data[[group_by]]))

# Create CellChat object using human ligand-receptor database
# (axolotl symbols are mapped to human via loc_map.tsv)
cellchat <- createCellChat(object = obj, group.by = group_by, assay = "SCT")
cellchat@DB <- CellChatDB.human

# Identify over-expressed signaling genes and interactions
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Infer communication probabilities at L-R level and filter low-confidence interactions
cellchat <- computeCommunProb(cellchat, type = "triMean", seed.use = 271)  # default: 1
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Aggregate to pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate communication network and compute centrality scores
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Aggregated interaction count circle plot
png(file.path(outdir, paste0(tp_id, "_interactions_count.png")),
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
png(file.path(outdir, paste0(tp_id, "_interactions_weight.png")),
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
png(file.path(outdir, paste0(tp_id, "_signaling_role.png")),
    width = 1600, height = 1200, res = 150)
print(netAnalysis_signalingRole_scatter(cellchat) +
      ggtitle(tp) +
      theme(aspect.ratio = 0.75))
dev.off()

# Pathway-by-cell-group heatmaps for outgoing and incoming signaling
png(file.path(outdir, paste0(tp_id, "_pathway_heatmap_outgoing.png")),
    width = 1200, height = 2400, res = 150)
print(netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",
                                        height = 30, font.size = 6))
dev.off()

png(file.path(outdir, paste0(tp_id, "_pathway_heatmap_incoming.png")),
    width = 1200, height = 2400, res = 150)
print(netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                        height = 30, font.size = 6))
dev.off()

# Save CellChat object
saveRDS(cellchat, file = file.path(outdir, paste0(tp_id, "_cellchat.rds")))

# Summary table (one-row TSV; merge step concatenates across timepoints)
write.table(
  data.frame(
    id                   = tp_id,
    timepoint            = tp,
    group_by             = group_by,
    n_cells              = length(cellchat@idents),
    n_groups             = length(levels(cellchat@idents)),
    n_signaling_genes    = length(cellchat@var.features$features),
    n_LR_pairs           = nrow(cellchat@LR$LRsig),
    n_significant_LR     = sum(cellchat@net$pval < 0.05, na.rm = TRUE),
    n_signaling_pathways = length(cellchat@netP$pathways)
  ),
  file      = file.path(outdir, paste0(tp_id, "_cellchat_summary.tsv")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)