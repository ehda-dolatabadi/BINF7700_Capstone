#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(SingleCellExperiment)
  library(scDblFinder) 
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id      <- args[1]
outdir  <- args[2]
input   <- args[3]

# Load mapped object
obj <- readRDS(input)

# Store initial counts
n_cells_before <- ncol(obj)
n_features_before <- nrow(obj)

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(obj)

# Run scDblFinder
sce <- scDblFinder(sce)

# Create vector of singlets
singlet_barcodes <- rownames(colData(sce))[sce$scDblFinder.class == "singlet"]

# Subset out doublets in the Seurat object
obj <- subset(obj, cells = singlet_barcodes)

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_DB_removed.rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells_before,
    n_cells_after = ncol(obj),
    n_features_before,
    n_features_after = nrow(obj)
  ),
  file = file.path(outdir, paste0(id, "_DB_removed_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# QC metrics
obj[["percent.mt"]]     <- PercentageFeatureSet(obj, pattern = "^(COX|ND|CYTB|ATP)")
obj[["percent.ribo"]]   <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
obj[["log10UMIsPerGene"]] <- log10(obj$nCount_RNA / obj$nFeature_RNA)

qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
df <- as_tibble(obj[[]], rownames="Cell.Barcode")

# Plots
pdf(file.path(outdir, paste0(id, "_DB_removed_plots.pdf")), width = 15, height = 9)

# QC plots
for (i in qc_metrics) {
  # scatter plots
  if (i != "nCount_RNA") {
    p <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = i) +
      theme(legend.title = element_blank())

    ggsave(file.path(outdir, paste0(id, "_", i, "_vs_nCount_RNA_DB_removed.png")), p,
           width = 15, height = 9, dpi = 300, bg = "white")
    print(p)
  }
}

dev.off()
