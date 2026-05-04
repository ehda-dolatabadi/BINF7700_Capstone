#!/usr/bin/env Rscript
# Script: merge_timepoints.R
# Purpose: Merge per-timepoint CellChat objects and produce cross-timepoint plots
# Description: Loads all per-timepoint CellChat RDS files written by the array job,
#              merges them with mergeCellChat, generates comparison plots
#              (total interaction count/strength, pathway information flow ranking),
#              and concatenates per-timepoint summary rows into a single TSV.
# Usage: Rscript merge_timepoints.R <id> <outdir>

# Load required libraries
suppressPackageStartupMessages({
  library(CellChat)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id     <- args[1]
outdir <- args[2]

# Canonical timepoint order from aggr_samples.csv
timepoint_order <- c("control", "3h", "24h", "72h", "7dpa", "14dpa", "22dpa", "33dpa")

# Discover per-timepoint outputs
cellchat_list <- list()
summary_rows  <- list()

for (tp in timepoint_order) {
  tp_id    <- paste0(id, "_", tp)
  rds_path <- file.path(outdir, tp_id, paste0(tp_id, "_cellchat.rds"))
  tsv_path <- file.path(outdir, tp_id, paste0(tp_id, "_cellchat_summary.tsv"))

  if (file.exists(rds_path)) {
    cellchat_list[[tp]] <- readRDS(rds_path)
    summary_rows[[tp]]  <- read.table(tsv_path, header = TRUE, sep = "\t",
                                      stringsAsFactors = FALSE)
  }
}

# Merge across timepoints and generate comparison plots
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
  # do.stat = FALSE because per-timepoint objects have different L-R pair sets
  # (significance testing requires matched pairs)
  png(file.path(outdir, paste0(id, "_rank_net.png")),
      width = 1200, height = 1800, res = 150)
  print(rankNet(cellchat_merged, mode = "comparison",
                stacked = TRUE, do.stat = FALSE))
  dev.off()

  # Save merged object
  saveRDS(cellchat_merged,
          file = file.path(outdir, paste0(id, "_cellchat_merged.rds")))
}

# Concatenated summary across timepoints
write.table(
  do.call(rbind, summary_rows),
  file      = file.path(outdir, paste0(id, "_cellchat_summary.tsv")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)