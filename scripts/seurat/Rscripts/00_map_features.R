#!/usr/bin/env Rscript
# Script: 00_map_features.R
# Purpose: Map gene IDs to gene symbols and create Seurat object from 10x data
# Description: Reads 10x CellRanger output, maps gene IDs to symbols using a reference TSV,
#              aggregates duplicate gene symbols, and creates a Seurat object
# Usage: Rscript 00_map_features.R <id> <outdir> <input> <tsv> <min_cells> <min_features>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(271)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id           <- args[1]
outdir       <- args[2]
input        <- args[3]
tsv          <- args[4]
min_cells    <- as.numeric(args[5])
min_features <- as.numeric(args[6])

# Load 10x data and create initial Seurat object
counts <- Read10X(input)
obj <- CreateSeuratObject(
  counts = counts,
  min.cells = min_cells,        # keep genes detected in ≥5 cells
  min.features = min_features,  # keep cells with ≥500 genes
  project = id
)

# Load gene ID to symbol mapping file and filter for genes in dataset
n_cells <- ncol(obj)
n_features_before <- nrow(obj)

loc_map <- read.table(tsv,
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE) %>% 
  filter(gene_id %in% rownames(obj))

n_genes_mapped <- nrow(loc_map)

# Extract count matrix and prepare for remapping
counts <- GetAssayData(obj, layer = "counts", assay = "RNA")

# Map gene IDs to gene symbols
new_names <- rownames(obj)
has_mapping <- rownames(obj) %in% loc_map$gene_id
new_names[has_mapping] <- loc_map$symbol[match(rownames(obj)[has_mapping], loc_map$gene_id)]
n_duplicates <- n_features_before - length(unique(new_names))

# Aggregate duplicate gene symbols by summing their counts
counts_agg <- rowsum(as.matrix(counts), group = new_names)

# Create new Seurat object with mapped gene symbols
obj <- CreateSeuratObject(
  counts = counts_agg,
  meta.data = obj@meta.data,
  project = id
)
n_features_after <- nrow(obj)

# Save mapped Seurat object
saveRDS(obj, file = file.path(dirname(outdir), paste0(id, ".rds")))

# Write summary statistics to file
write.table(
  data.frame(
    id,
    min_cells,
    min_features,
    n_cells,
    n_features_before,
    n_features_after,
    n_genes_mapped,
    n_duplicates
  ),
  file = file.path(outdir, paste0(id, "_mapping_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
