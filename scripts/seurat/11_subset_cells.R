#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(future)
})

set.seed(777)

# Set up parallelization
options(future.globals.maxSize = 16000 * 1024^2)  # 16 GB
plan("multicore", workers = 16)

# Close any open graphics devices
while (!is.null(dev.list())) dev.off()

# Args
args <- commandArgs(trailingOnly = TRUE)
id	<- args[1]
outdir	<- args[2]
input	<- args[3]

# Parameters
epi_thr <- 0.5
ery_thr <- 3.0

# Load clustered object
obj <- readRDS(input)

# Use SCT assay
DefaultAssay(obj) <- "SCT"

n_cells_before = ncol(obj)
n_features_before = nrow(obj)

# Define marker sets
epithelial_markers <- c("KRT8", "KRT18", "KRT1", "EPCAM", "CLDN7", 
                        "CDH1", "DSP", "DSG2", "EPPK1", "S100P")

erythrocyte_markers <- c("HBA1", "HBB", "HBA2", "HEMGN", 
                         "ALAS2", "PRDX2", "ANK1", "HBBP1")

epithelial_markers_available <- epithelial_markers[epithelial_markers %in% rownames(obj)]
erythrocyte_markers_available <- erythrocyte_markers[erythrocyte_markers %in% rownames(obj)]

# Add module score for epithelial cells
obj <- AddModuleScore(
  obj,
  features = list(epithelial_markers_available),
  name = "Epithelial_score"
)
n_cells_epithelial = sum(obj$Epithelial_score1 >= epi_thr)

# Add module score for erythrocyte cells
obj <- AddModuleScore(
  obj,
  features = list(erythrocyte_markers_available),
  name = "Erythrocyte_score"
)
n_cells_erythrocyte = sum(obj$Erythrocyte_score1 >= ery_thr)

# Plot VlnPlot
png(file.path(outdir, paste0(id, "_00_violin_score.png")), width = 1600, height = 1200)
VlnPlot(obj, features = c("Epithelial_score1", "Erythrocyte_score1"), pt.size = 0, group.by = "orig.ident") +
  labs(x = NULL)
dev.off()

# Visualize distribution
png(file.path(outdir, paste0(id, "_00_epithelial_score_distribution.png")), width = 1600, height = 1200)
hist(obj$Epithelial_score1, breaks = 50, main = "Epithelial Score Distribution", xlab = "Epithelial Score")
abline(v = epi_thr, col = "red", lwd = 2, lty = 2)
dev.off()

png(file.path(outdir, paste0(id, "_00_erythrocyte_score_distribution.png")), width = 1600, height = 1200)
hist(obj$Erythrocyte_score1, breaks = 50, main = "Erythrocyte Score Distribution", xlab = "Erythrocyte Score")
abline(v = ery_thr, col = "red", lwd = 2, lty = 2)
dev.off()

# Filtering
obj <- subset(
  obj,
  subset =
    Epithelial_score1 < epi_thr &
    Erythrocyte_score1 < ery_thr
)

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_00_subset.rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells_before,
    n_cells_after = ncol(obj),
    n_features_before,
    n_features_after = nrow(obj),
    n_epithelial_markers_total = length(epithelial_markers),
    n_epithelial_markers_available = length(epithelial_markers_available),
    n_erythrocyte_markers_total = length(erythrocyte_markers),
    n_erythrocyte_markers_available = length(erythrocyte_markers_available),
    epi_thr,
    ery_thr,
    n_cells_epithelial,
    pct_cells_epithelial = n_cells_epithelial / n_cells_before * 100,
    n_cells_erythrocyte,
    pct_cells_erythrocyte = n_cells_erythrocyte / n_cells_before * 100,
    pct_cells_kept = ncol(obj) / n_cells_before * 100
  ),
  file = file.path(outdir, paste0(id, "_00_subsetting_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
