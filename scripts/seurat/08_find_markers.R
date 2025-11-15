#!/usr/bin/env Rscript
# Purpose: Find marker genes for each Seurat cluster
# Usage: Rscript 8_find_markers.R <id> <output> <input> <cores>

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]
cores	<- args[4]

# Set up parallelization
options(future.globals.maxSize = 16000 * 1024^2)  # 16 GB
plan("multicore", workers = as.numeric(cores))

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
)

# Order markers
markers <- markers[order(markers$cluster, markers$p_val_adj, -markers$avg_log2FC), ]

# Save markers
outfile <- file.path(outdir, paste0(id, "_markers.tsv"))
write.table(markers, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)

# Filter
markers <- markers %>%
  mutate(direction = ifelse(avg_log2FC >= 0, "up", "down")) %>%
  filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) >= 1,
    abs(pct.1 - pct.2) >= 0.2
  )

# Rank and pick top 10 per cluster and direction
top <- markers %>%
  group_by(cluster, direction) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)), desc(abs(pct.1 - pct.2))) %>%
  slice_head(n = 10) %>%
  ungroup()

# Save filtered markers
outfile <- file.path(outdir, paste0(id, "_markers_filtered.tsv"))
write.table(top, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
