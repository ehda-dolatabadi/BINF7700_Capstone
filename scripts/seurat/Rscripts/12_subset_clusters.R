#!/usr/bin/env Rscript
# Script: 13_subset_clusters.R
# Purpose: Extract cells by cluster identity and SingleR label
# Usage: Rscript 13_subset_clusters.R <id> <outdir> <input> <mode> <value>
#
#   mode : keep or remove
#   value  : cluster ID(s) and SingleR label(s) (semicolon-separated)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(271)

args <- commandArgs(trailingOnly = TRUE)
id     <- args[1]
outdir <- args[2]
input  <- args[3]
mode   <- args[4]
idents  <- unlist(strsplit(args[5], ";"))
labels  <- unlist(strsplit(args[6], ";"))

if (!mode %in% c("keep", "remove")) {
  stop("mode must be 'keep' or 'remove'")
}

options(future.seed = TRUE)
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 56)

obj <- readRDS(input)
n_cells_before    <- ncol(obj)
n_features_before <- nrow(obj[["SCT"]])

if (mode == "keep") {
  obj <- subset(obj, seurat_clusters %in% idents | singler_label %in% labels)
} else {
  obj <- subset(obj, !(seurat_clusters %in% idents & singler_label %in% labels))
}

saveRDS(obj, file = file.path(dirname(outdir), paste0(id, "_integrated.rds")))

write.table(
  data.frame(
    id,
    mode,
    idents = paste0(idents, collapse=","),
	labels = paste0(labels, collapse=","),
    n_cells_before,
    n_cells_after     = ncol(obj),
    n_features_before,
    n_features_after  = nrow(obj[["SCT"]])
  ),
  file      = file.path(outdir, paste0(id, "_subsetting_summary.tsv")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)