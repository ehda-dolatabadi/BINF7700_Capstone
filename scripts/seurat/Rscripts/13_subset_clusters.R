#!/usr/bin/env Rscript
# Script: 13_subset_clusters.R
# Purpose: Extract cells by cluster identity or SingleR label
# Description: Subsets the Seurat object to retain only cells matching the
#              specified cluster idents or singler_label for downstream analysis
# Usage: Rscript 13_subset_clusters.R <id> <outdir> <input> <method> <value>
#
#   method : "idents"  — subset by Seurat cluster identity (seurat_clusters)
#            "singler" — subset by SingleR label (singler_label metadata column)
#   value  : cluster ID(s) or SingleR label to keep (semicolon-separated for multiple)

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
method <- args[4]
value  <- args[5]

if (!method %in% c("idents", "singler")) {
  stop("method must be 'idents' or 'singler'")
}

options(future.seed = TRUE)
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 56)

obj <- readRDS(input)
n_cells_before    <- ncol(obj)
n_features_before <- nrow(obj[["SCT"]])

value_vec <- unlist(strsplit(value, ";"))

if (method == "idents") {
  obj <- subset(obj, idents = value_vec)
} else {
  obj <- JoinLayers(obj, assay = "RNA")
  obj <- subset(obj, singler_label %in% value_vec)
}

saveRDS(obj, file = file.path(dirname(outdir), paste0(id, "_integrated.rds")))

write.table(
  data.frame(
    id,
    method,
    value,
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