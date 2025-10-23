#!/usr/bin/env Rscript
# Purpose: Run QC on a single 10x sample

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript qc.R <id> <input> <outdir>")
}

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
seurat <- CreateSeuratObject(
  counts = counts,
  min.cells = min_cells,	# keep genes detected in ≥5 cells
  min.features = min_features,	# keep cells with ≥500 genes
  project = id
)

# QC metrics
seurat[["log10UMIsPerGene"]]	<- log10(seurat$nCount_RNA / seurat$nFeature_RNA)
seurat[["percent.mt"]]		<- PercentageFeatureSet(seurat, pattern = "^(COX|ND|CYTB|ATP)")
seurat[["percent.ribo"]]	<- PercentageFeatureSet(seurat, pattern = "^RPS|^RPL")
seurat[["percent.rrna"]]	<- PercentageFeatureSet(seurat, pattern = "^RRN")
seurat[["percent.trna"]]	<- PercentageFeatureSet(seurat, pattern = "^TRNA")
seurat[["percent.loc"]]		<- PercentageFeatureSet(seurat, pattern = "^LOC")

qc_metrics <- c("log10UMIsPerGene","percent.mt","percent.ribo","percent.rrna","percent.trna","percent.loc")

# Plots
png(file.path(outdir, paste0(id, "_feature_vs_count.png")), width=1200, height=720)
print(FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

for (i in qc_metrics) {
  png(file.path(outdir, paste0(id, "_", i, ".png")), width=1200, height=720)
  print(VlnPlot(seurat, features = i) + NoLegend() + labs(x = NULL))
  dev.off()
}

# Filters
n_before <- ncol(seurat)
seurat_filter <- subset(
  seurat,
  subset =
    nFeature_RNA <= max_features &
    nCount_RNA   >= min_counts   &
    nCount_RNA   <= max_counts   &
    percent.mt   <= max_mt       &
    percent.ribo <= max_ribo
)
n_after <- ncol(seurat_filter)

# Save
saveRDS(seurat_filter, file = file.path(outdir, paste0(id, "_filtered.rds")))

write.table(
  data.frame(id, n_cells_before=n_before, n_cells_after=n_after,
             min_cells, min_features, max_features,
	     min_counts, max_counts,
             max_mt, max_ribo),
  file = file.path(outdir, paste0(id, "_qc_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(sessionInfo())
