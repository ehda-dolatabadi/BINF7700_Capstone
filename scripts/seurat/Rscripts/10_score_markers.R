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
  name = cell_name
)
feature_name = paste0(cell_name, "1")

# Plot VlnPlot
png(file.path(outdir, paste0(id, "_all_timepoints_violin_score.png")), width = 1600, height = 1200, res = 150)
print(VlnPlot(obj, features = feature_name, group.by = "orig.ident", pt.size = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = paste0(cell_name, " by Timepoint"),
    x = "Timepoint",
    y = "Module Score"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  ))
dev.off()

png(file.path(outdir, paste0(id, "_all_clusters_violin_score.png")), width = 1600, height = 1200, res = 150)
print(VlnPlot(obj, features = feature_name, group.by = "seurat_clusters", pt.size = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = paste0(cell_name, " by Cluster"),
    x = "Cluster",
    y = "Module Score"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  ))
dev.off()

# Plot individual markers
for (marker in markers_available){
  png(file.path(outdir, paste0(id, "_", marker, "_timepoints_violin_score.png")), width = 1600, height = 1200, res = 150)
  print(VlnPlot(obj, features = marker, group.by = "orig.ident", pt.size = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(
      title = paste0(marker, " Expression by Timepoint"),
      x = "Timepoint",
      y = "Expression Level"
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11)
    ))
  dev.off()
}

for (marker in markers_available){
  png(file.path(outdir, paste0(id, "_", marker, "_clusters_violin_score.png")), width = 1600, height = 1200, res = 150)
  print(VlnPlot(obj, features = marker, group.by = "seurat_clusters", pt.size = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(
      title = paste0(marker, " Expression by Cluster"),
      x = "Cluster",
      y = "Expression Level"
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11)
    ))
  dev.off()
}

# Visualize distribution
png(file.path(outdir, paste0(id, "_all_score_distribution.png")), width = 1600, height = 1200, res = 150)
hist(obj@meta.data[[feature_name]], breaks = 50,
     main = paste0(cell_name, " Distribution"),
     xlab = "Module Score",
     ylab = "Frequency (Number of Cells)",
     col = "lightblue",
     border = "black",
     cex.main = 1.5,
     cex.lab = 1.3,
     cex.axis = 1.1)
abline(v = 0, col = "gray", lty = 2, lwd = 2)
legend("topright", legend = "Zero reference", lty = 2, col = "gray", lwd = 2, cex = 1.1)
dev.off()

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_scored.rds")))

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
