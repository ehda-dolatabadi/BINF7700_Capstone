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
png(file.path(outdir, paste0(id, "_timepoints_violin_score.png")), width = 1600, height = 1200, res = 150)
VlnPlot(obj, features = feature_name, group.by = "orig.ident", pt.size = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = NULL)
dev.off()

# Plot individual markers
for (marker in markers_available){
  png(file.path(outdir, paste0(id, marker, "_timepoints_violin_score.png")), width = 1600, height = 1200, res = 150)
  print(VlnPlot(obj, features = marker, group.by = "orig.ident", pt.size = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(x = NULL))
  dev.off()
}

if ("seurat_clusters" %in% colnames(obj@meta.data)) {
  png(file.path(outdir, paste0(id, "_clusters_violin_score.png")), width = 1600, height = 1200, res = 150)
  VlnPlot(obj, features = feature_name, group.by = "seurat_clusters", pt.size = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(x = NULL)
  dev.off()

  for (marker in markers_available){
    png(file.path(outdir, paste0(id, marker, "_clusters_violin_score.png")), width = 1600, height = 1200, res = 150)
    print(VlnPlot(obj, features = marker, group.by = "seurat_clusters", pt.size = 0) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(x = NULL))
    dev.off()
  }
}

# Visualize distribution
png(file.path(outdir, paste0(id, "_score_distribution.png")), width = 1600, height = 1200, res = 150)
hist(obj@meta.data[[feature_name]], breaks = 50, main = "Score Distribution", xlab = "Score")
abline(v = 0, col = "red", lty = 2, lwd = 2)
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
