#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]
tsv	<- args[4]

# Parameters
min_cells       <- 5
min_features    <- 500

# Load and create Seurat object
counts <- Read10X(input)
obj <- CreateSeuratObject(
  counts = counts,
  min.cells = min_cells,        # keep genes detected in ≥5 cells
  min.features = min_features,  # keep cells with ≥500 genes
  project = id
)

# Load mapping file
loc_map <- read.table(tsv,
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE) %>% 
  filter(gene_id %in% rownames(obj))

# Get counts
counts <- GetAssayData(obj, slot = "counts", assay = "RNA")

# Map ALL rownames
new_names <- rownames(obj)

# Apply mapping
has_mapping <- rownames(obj) %in% loc_map$gene_id
new_names[has_mapping] <- loc_map$symbol[match(rownames(obj)[has_mapping], loc_map$gene_id)]

# Aggregate by summing duplicates
counts_agg <- rowsum(as.matrix(counts), group = new_names)

# Create new object
obj <- CreateSeuratObject(
  counts = counts_agg,
  meta.data = obj@meta.data
)

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_mapped.rds")))
