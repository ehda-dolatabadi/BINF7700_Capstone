#!/usr/bin/env Rscript
# Purpose: Run filtering on a 10x sample
# Usage: Rscript 2_filter.R <id> <output> <input>

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id           <- args[1]
outdir       <- args[2]
input        <- args[3]
min_counts   <- as.numeric(args[4])
max_counts   <- as.numeric(args[5])
min_features <- as.numeric(args[6])
max_features <- as.numeric(args[7])
max_mt       <- as.numeric(args[8])
max_ribo     <- as.numeric(args[9])

# Load mapped object
obj <- readRDS(input)
n_cells_before <- ncol(obj)
n_features_before <- nrow(obj)

# QC metrics
obj[["percent.mt"]]     <- PercentageFeatureSet(obj, pattern = "^(COX|ND|CYTB|ATP)")
obj[["percent.ribo"]]   <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")

qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt","percent.ribo")
df <- as_tibble(obj[[]], rownames="Cell.Barcode")

# Filter
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

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_filtered.rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells_before,
    n_cells_after = ncol(obj),
    n_features_before,
    n_features_after = nrow(obj),
    min_counts,
    max_counts,
    min_features,
    max_features,
    max_mt,
    max_ribo
  ),
  file = file.path(outdir, paste0(id, "_filtering_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Plots
pdf(file.path(outdir, paste0(id, "_filtered_plots.pdf")), width = 15, height = 9)
for (i in qc_metrics) {

  # scatter plots
  if (i != "nCount_RNA") {
    p <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = i) +
      theme(legend.title = element_blank())

    ggsave(file.path(outdir, paste0(id, "_", i, "_vs_nCount_RNA_filtered.png")), p,
           width = 15, height = 9, dpi = 300, bg = "white")
    print(p)
  }
}

dev.off()
