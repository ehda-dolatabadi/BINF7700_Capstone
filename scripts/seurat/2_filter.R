#!/usr/bin/env Rscript
# Purpose: Run filtering on a 10x sample
# Usage: Rscript 2_filter.R <id> <input> <outdir>

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)

id	<- args[1]
input	<- args[2]
outdir	<- args[3]

# Thresholds
min_cells	<- 5
min_features	<- 500
max_features	<- 7500
min_counts	<- 1000
max_counts	<- 40000
max_mt		<- 25
max_ribo	<- 40

# Load and create Seurat object
counts <- Read10X(input)
obj <- CreateSeuratObject(
  counts = counts,
  min.cells = min_cells,	# keep genes detected in ≥5 cells
  min.features = min_features,	# keep cells with ≥500 genes
  project = id
)

# Filters
n_before <- ncol(obj)
obj <- subset(
  obj,
  subset =
    nFeature_RNA <= max_features &
    nCount_RNA   >= min_counts   &
    nCount_RNA   <= max_counts   &
    percent.mt   <= max_mt       &
    percent.ribo <= max_ribo
)
n_after <- ncol(obj)

# Save
saveRDS(obj, file = file.path(outdir, paste0(id, "_filtered.rds")))

write.table(
  data.frame(id, n_cells_before=n_before, n_cells_after=n_after,
             min_cells, min_features, max_features,
	     min_counts, max_counts,
             max_mt, max_ribo),
  file = file.path(outdir, paste0(id, "_filtered.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(sessionInfo())
