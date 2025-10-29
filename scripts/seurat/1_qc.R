#!/usr/bin/env Rscript
# Purpose: Run QC on a 10x sample
# Usage: Rscript 1_qc.R <id> <output> <input>

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id	<- args[1]
outdir	<- args[2]
input	<- args[3]

# Parameters
min_cells	<- 5
min_features	<- 500

# Load and create Seurat object
counts <- Read10X(input)
obj <- CreateSeuratObject(
  counts = counts,
  min.cells = min_cells,	# keep genes detected in ≥5 cells
  min.features = min_features,	# keep cells with ≥500 genes
  project = id
)

# QC metrics
obj[["log10UMIsPerGene"]] <- log10(obj$nCount_RNA / obj$nFeature_RNA)

obj[["percent.mt"]]	<- PercentageFeatureSet(obj, pattern = "^(COX|ND|CYTB|ATP)")
obj[["percent.ribo"]]	<- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
obj[["percent.rrna"]]	<- PercentageFeatureSet(obj, pattern = "^RRN")
obj[["percent.trna"]]	<- PercentageFeatureSet(obj, pattern = "^TRNA")
obj[["percent.loc"]]	<- PercentageFeatureSet(obj, pattern = "^LOC")

qc_metrics <- c("log10UMIsPerGene","percent.mt","percent.ribo","percent.rrna","percent.trna","percent.loc")

# Plots
png(file.path(outdir, paste0(id, "_feature_vs_count.png")), width=1200, height=720)
print(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

for (i in qc_metrics) {
  png(file.path(outdir, paste0(id, "_", i, ".png")), width=1200, height=720)
  print(VlnPlot(obj, features = i) + NoLegend() + labs(x = NULL))
  dev.off()
}
