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

# Load clustered object
obj <- readRDS(input)

# Use SCT assay
DefaultAssay(obj) <- "SCT"

# Core Schwann cell identity markers (highest specificity)
schwann_markers <- c("SOX10", "S100", "S100B", "NGFR", "p75NTR", "MPZ", "MBP", "PMP22", "PLP1", "PRX", 
		     "NCAM", "NCAM1", "L1CAM", "SCN7A", "SOX2", "GAP43", "EGR2", "Krox20", "POU3F1", "OCT6")

schwann_markers_available <- schwann_markers[schwann_markers %in% rownames(obj)]

# Add module score (makes schwann1 feature)
obj <- AddModuleScore(
  obj,
  features = list(schwann_markers_available),
  name = "Schwann"
)

# Plot FeaturePlot
png(file.path(outdir, paste0(id, "_feature_score.png")), width = 1600, height = 1200)
FeaturePlot(obj, features = "Schwann1")
dev.off()

# Plot VlnPlot
png(file.path(outdir, paste0(id, "_violin_score.png")), width = 1600, height = 1200)
VlnPlot(obj, features = "Schwann1", group.by = "seurat_clusters", pt.size = 0) +
  labs(x = NULL)
dev.off()

# Check individual marker expression
png(file.path(outdir, paste0(id, "_individual_markers_violin.png")), width = 1600, height = 1200)
VlnPlot(obj, features = schwann_markers_available, ncol = 2, pt.size = 0) +
  labs(x = NULL)
dev.off()

# Visualize distribution
png(file.path(outdir, paste0(id, "_schwann_score_distribution.png")), width = 1600, height = 1200)
hist(obj$Schwann1, breaks = 50, main = "Schwann Score Distribution", xlab = "Schwann Score")
abline(v = quantile(obj$Schwann1, 0.9), col = "red", lty = 2)  # Top 10%
abline(v = quantile(obj$Schwann1, 0.95), col = "blue", lty = 2) # Top 5%
legend("topright", c("90th percentile", "95th percentile"), col = c("red", "blue"), lty = 2)
dev.off()

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_scored.rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells = ncol(obj),
    n_features = nrow(obj),
    n_schwann_markers_total = length(schwann_markers),
    n_schwann_markers_available = length(schwann_markers_available),
    schwann_markers_available = paste(schwann_markers_available, collapse = ",")
  ),
  file = file.path(outdir, paste0(id, "_schwann_score_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
