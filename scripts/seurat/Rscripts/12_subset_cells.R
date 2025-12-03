#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(future)
})

set.seed(777)

# Set up parallelization
options(future.globals.maxSize = 64000 * 1024^2)  # 64 GB
plan("multicore", workers = 56)

# Close any open graphics devices
while (!is.null(dev.list())) dev.off()

# Args
args <- commandArgs(trailingOnly = TRUE)
id			<- args[1]
outdir			<- args[2]
input			<- args[3]
mode			<- args[4]
cell_types_str		<- args[5]
cell_markers_str	<- args[6]
cell_thresholds_str	<- args[7]

# Parse arguments
cell_types <- unlist(strsplit(cell_types_str, ";"))
cell_markers_list <- strsplit(unlist(strsplit(cell_markers_str, ";")), ",")
cell_thresholds <- as.numeric(unlist(strsplit(cell_thresholds_str, ";")))

# Validate inputs
if (length(cell_types) != length(cell_markers_list) || length(cell_types) != length(cell_thresholds)) {
  stop("Number of cell types, marker lists, and thresholds must match")
}

# Load clustered object
obj <- readRDS(input)

# Use SCT assay
DefaultAssay(obj) <- "SCT"

n_cells_before = ncol(obj)
n_features_before = nrow(obj)

# Initialize summary data
summary_data <- data.frame(
  id = id,
  n_cells_before = n_cells_before,
  n_features_before = n_features_before
)

# Process each cell type
score_features <- c()

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  markers <- cell_markers_list[[i]]
  threshold <- cell_thresholds[i]

  # Find available markers
  markers_available <- markers[markers %in% rownames(obj)]

  # Add module score
  score_name <- paste0(cell_type, "_score")
  obj <- AddModuleScore(
    obj,
    features = list(markers_available),
    name = score_name
  )

  # Note: AddModuleScore adds "1" to the name
  score_name <- paste0(score_name, "1")
  score_features <- c(score_features, score_name)

  # Count cells above threshold
  n_cells_above <- sum(obj[[score_name]] >= threshold)

  # Distribution plot
  png(file.path(outdir, paste0(id, "_", cell_type, "_score_distribution.png")),
      width = 1600, height = 1200, res=150)
  hist(obj[[score_name]][,1], breaks = 50,
       main = paste(cell_type, "Score Distribution"),
       xlab = paste(cell_type, "Score"))
  abline(v = c(0, threshold), col = c("gray", "red"), lwd = 2, lty = 2)
  dev.off()
  
  # Violin plot
  png(file.path(outdir, paste0(id, "_", cell_type, "_violin_score.png")), width = 1600, height = 1200, res=150)
  print(VlnPlot(obj, features = score_name, pt.size = 0, group.by = "orig.ident") +
    geom_hline(yintercept = c(0, threshold), linetype = "dashed", color = c("gray", "red")) +
    labs(x = NULL))
  dev.off()

  # Add to summary
  summary_data[[paste0("n_", cell_type, "_markers_total")]] <- length(markers)
  summary_data[[paste0("n_", cell_type, "_markers_available")]] <- length(markers_available)
  summary_data[[paste0(cell_type, "_threshold")]] <- threshold
  summary_data[[paste0("n_cells_", cell_type)]] <- n_cells_above
  summary_data[[paste0("pct_cells_", cell_type)]] <- n_cells_above / n_cells_before * 100
}

# To remove
if (mode == "remove") {
  # Get cells that pass the filter
  cells_to_keep <- rep(TRUE, ncol(obj))
  for (i in 1:length(cell_types)) {
    score_col <- paste0(cell_types[i], "_score1")
    threshold <- cell_thresholds[i]
    cells_to_keep <- cells_to_keep & (obj[[score_col]][,1] < threshold)
  }
  # Subset to keep only cells that pass all filters
  obj <- obj[, cells_to_keep]
}

# To enrich
if (mode == "keep") {
  # Get cells that pass the filter
  cells_to_remove <- rep(TRUE, ncol(obj))
  for (i in 1:length(cell_types)) {
    score_col <- paste0(cell_types[i], "_score1")
    threshold <- cell_thresholds[i]
    cells_to_remove <- cells_to_remove & (obj[[score_col]][,1] < threshold)
  }
  # Subset to keep only cells that pass all filters
  obj <- obj[, !cells_to_remove]
}

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_subset.rds")))

# Complete summary
summary_data$n_cells_after <- ncol(obj)
summary_data$n_features_after <- nrow(obj)
summary_data$pct_cells_kept <- ncol(obj) / n_cells_before * 100

# Save summary
write.table(
  summary_data,
  file = file.path(outdir, paste0(id, "_subsetting_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
