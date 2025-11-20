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

# Set up parallelization
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 64)

# Parameters
significance <- 0.05
regulation <- 1
enrichment <- 0.2
top <- 30

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
outfile <- file.path(outdir, paste0(id, "_08_markers.tsv"))
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
