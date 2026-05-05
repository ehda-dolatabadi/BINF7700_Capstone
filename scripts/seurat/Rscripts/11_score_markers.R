#!/usr/bin/env Rscript
# Script: 11_score_markers.R
# Purpose: Calculate module scores for cell type markers and visualize
# Description: Computes module scores based on predefined marker genes to identify
#              cell populations, generates violin plots and expression plots
# Usage: Rscript 11_score_markers.R <id> <outdir> <input> <cell_name> <markers>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(271)

# Close any open graphics devices
while (!is.null(dev.list())) dev.off()

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
outdir <- args[2]
input <- args[3]
cell_name <- args[4]
markers <- unlist(strsplit(args[5], ","))

# Load clustered object
obj <- readRDS(input)

# Use SCT assay
DefaultAssay(obj) <- "SCT"

# Detect markers
markers_available <- markers[markers %in% rownames(obj)]
markers_missing <- markers[!markers %in% rownames(obj)]

if (length(markers_available) == 0) {
  cat("None of the markers found in dataset.")
  quit(save = "no", status = 0)
}

# Add module score
cell_name <- paste0(cell_name, "_score")
obj <- AddModuleScore(
  obj,
  features = list(markers_available),
  name = cell_name,
  seed = 271  # default: 1
)
feature_name = paste0(cell_name, "1")

# Plot VlnPlot
png(
  file.path(outdir, paste0(id, "_all_timepoints_violin_score.png")),
  width = 1600, height = 1200, res = 150
)
print(VlnPlot(obj, features = feature_name, group.by = "orig.ident", pt.size = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = cell_name,
    x = "Timepoint",
    y = "Module Score"
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  ))
dev.off()

png(
  file.path(outdir, paste0(id, "_all_clusters_violin_score.png")),
  width = 1600, height = 1200, res = 150
)
print(VlnPlot(obj, features = feature_name, group.by = "seurat_clusters", pt.size = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = cell_name,
    x = "Cluster",
    y = "Module Score"
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  ))
dev.off()

# Plot individual markers
for (marker in markers_available){
  png(
    file.path(outdir, paste0(id, "_", marker, "_timepoints_violin_score.png")),
    width = 1600, height = 1200, res = 150
  )
  print(VlnPlot(obj, features = marker, group.by = "orig.ident", pt.size = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(
      title = marker,
      x = "Timepoint",
      y = "Expression Level"
    ) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14)
    ))
  dev.off()
}

for (marker in markers_available){
  png(
    file.path(outdir, paste0(id, "_", marker, "_clusters_violin_score.png")),
    width = 1600, height = 1200, res = 150
  )
  print(VlnPlot(obj, features = marker, group.by = "seurat_clusters", pt.size = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(
      title = marker,
      x = "Cluster",
      y = "Expression Level"
    ) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14)
    ))
  dev.off()
}

# Visualize distribution
png(
  file.path(outdir, paste0(id, "_all_score_distribution.png")),
  width = 1600, height = 1200, res = 150
)
hist(obj@meta.data[[feature_name]], breaks = 50,
     main = paste0(cell_name, " Distribution"),
     xlab = "Module Score",
     ylab = "Frequency (Number of Cells)",
     col = "lightblue",
     border = "black",
     cex.main = 2.0,
     cex.lab = 1.8,
     cex.axis = 1.5)
abline(v = 0, col = "gray", lty = 2, lwd = 2)
legend("topright", legend = "Zero reference", lty = 2, col = "gray", lwd = 2, cex = 1.1)
dev.off()

# Summary table
write.table(
  data.frame(
    id,
    n_markers_total = length(markers),
    n_markers_available = length(markers_available),
    markers_available = paste(markers_available, collapse = ","),
    markers_missing = paste(markers_missing, collapse = ",")
  ),
  file = file.path(outdir, paste0(id, "_score_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
