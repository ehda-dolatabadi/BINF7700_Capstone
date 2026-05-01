#!/usr/bin/env Rscript
# Script: 09_find_markers.R
# Purpose: Identify differentially expressed marker genes for each cluster
# Description: Finds cluster-specific markers using Wilcoxon test, filters by significance
#              and expression thresholds, and exports top markers per cluster
# Usage: Rscript 09_find_markers.R <id> <output> <input>
#        <significance> <regulation> <enrichment> <top>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(271)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id           <- args[1]
outdir       <- args[2]
input        <- args[3]
significance <- as.numeric(args[4])
regulation   <- as.numeric(args[5])
enrichment   <- as.numeric(args[6])
top          <- as.numeric(args[7])

# Set up parallelization
options(future.seed = TRUE)
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 56)

# Load clustered object
obj <- readRDS(input)

# Use SCT assay and prepare it
DefaultAssay(obj) <- "SCT"
obj <- PrepSCTFindMarkers(obj)

# Find markers
markers <- FindAllMarkers(
  obj,
  max.cells.per.ident = 500,
  only.pos = FALSE,
  random.seed = 271  # default: 1; used by set.seed() internally before downsampling
)

# Order markers
markers <- markers[order(markers$cluster, markers$p_val_adj, -markers$avg_log2FC), ]

# Save markers
outfile <- file.path(outdir, paste0(id, "_markers.tsv"))
write.table(markers, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)

# Filter
filtered <- markers %>%
  mutate(direction = ifelse(avg_log2FC >= 0, "up", "down")) %>%
  filter(
    p_val_adj < significance,
    abs(avg_log2FC) >= regulation,
    abs(pct.1 - pct.2) >= enrichment
  )

# Rank and pick top rankings per cluster and direction
filtered_top <- filtered %>%
  group_by(cluster, direction) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)), desc(abs(pct.1 - pct.2))) %>%
  slice_head(n = top) %>%
  ungroup()

# Save filtered markers
outfile <- file.path(outdir, paste0(id, "_markers_filtered.tsv"))
write.table(filtered_top, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
