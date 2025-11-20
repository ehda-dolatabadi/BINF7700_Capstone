#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(777)

# Close any open graphics devices
while (!is.null(dev.list())) dev.off()

# Args
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
outdir <- args[2]
input <- args[3]

# Parameters
cell_name <- "Schwann"
markers <- c("SOX10", "S100", "S100B", "NGFR", "p75NTR", "MPZ", "MBP", "PMP22", "PLP1", "PRX",
             "NCAM", "NCAM1", "L1CAM", "SCN7A", "SOX2", "GAP43", "EGR2", "Krox20", "POU3F1", "OCT6")
group_by <- "seurat_clusters"

# Load clustered object
obj <- readRDS(input)

# Use SCT assay
DefaultAssay(obj) <- "SCT"

# Detect markers
markers_available <- markers[markers %in% rownames(obj)]

# Add module score
obj <- AddModuleScore(
  obj,
  features = list(markers_available),
  name = cell_name
)
feature_name = paste0(cell_name, "1")

# Plot FeaturePlot
png(file.path(outdir, paste0(id, "_09_feature_score.png")), width = 1600, height = 1200)
FeaturePlot(obj, features = feature_name)
dev.off()

# Plot VlnPlot
png(file.path(outdir, paste0(id, "_09_violin_score.png")), width = 1600, height = 1200)
VlnPlot(obj, features = feature_name, group.by = group_by, pt.size = 0) +
  labs(x = NULL)
dev.off()

# Check individual marker expression
png(file.path(outdir, paste0(id, "_09_individual_markers_violin.png")), width = 1600, height = 1200)
VlnPlot(obj, features = markers_available, ncol = 2, pt.size = 0) +
  labs(x = NULL)
dev.off()

# Visualize distribution
png(file.path(outdir, paste0(id, "_09_score_distribution.png")), width = 1600, height = 1200)
hist(obj[[feature_name]], breaks = 50, main = "Score Distribution", xlab = "Score")
dev.off()

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_09_scored.rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells = ncol(obj),
    n_features = nrow(obj),
    n_markers_total = length(markers),
    n_markers_available = length(markers_available),
    markers_available = paste(markers_available, collapse = ",")
  ),
  file = file.path(outdir, paste0(id, "_09_score_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
