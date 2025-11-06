#!/usr/bin/env Rscript
# Purpose: Run filtering on a 10x sample
# Usage: Rscript 2_filter.R <id> <output> <input> <params>

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]

# Parameters
min_cells    <- 5
min_counts   <- as.numeric(args[4])
max_counts   <- as.numeric(args[5])
min_features <- as.numeric(args[6])
max_features <- as.numeric(args[7])
max_mt       <- as.numeric(args[8])
max_ribo     <- as.numeric(args[9])

# Load and create Seurat object
counts <- Read10X(input)
obj <- CreateSeuratObject(
  counts = counts,
  min.cells = min_cells,	# keep genes detected in ≥5 cells
  min.features = min_features,	# keep cells with ≥500 genes
  project = id
)

# QC metrics
obj[["percent.mt"]]	<- PercentageFeatureSet(obj, pattern = "^(COX|ND|CYTB|ATP)")
obj[["percent.ribo"]]	<- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")

# Filter
n_before <- ncol(obj)
obj <- subset(
  obj,
  subset =
    nCount_RNA   >= min_counts   &
    nCount_RNA   <= max_counts   &
    nFeature_RNA >= min_features &
    nFeature_RNA <= max_features &
    percent.mt   <= max_mt       &
    percent.ribo <= max_ribo
)
n_after <- ncol(obj)

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_filtered.rds")))

# Summary
write.table(
  data.frame(id, n_cells_before=n_before, n_cells_after=n_after,
             min_cells, min_features, max_features,
	     min_counts, max_counts,
             max_mt, max_ribo),
  file = file.path(outdir, paste0(id, "_filtered.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
