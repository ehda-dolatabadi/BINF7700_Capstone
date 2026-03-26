#!/usr/bin/env Rscript
# Script: 03_filter.R
# Purpose: Filter cells based on QC thresholds
# Description: Applies filtering thresholds for UMI counts, feature counts, mitochondrial
#              and ribosomal content, and generates QC plots with threshold lines
# Usage: Rscript 03_filter.R <id> <output> <input> <min_counts> <max_counts> <min_features> <max_features> <max_mito> <max_ribo>

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(271)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
id           <- args[1]
outdir       <- args[2]
input        <- args[3]
min_counts   <- as.numeric(args[4])
max_counts   <- as.numeric(args[5])
min_features <- as.numeric(args[6])
max_features <- as.numeric(args[7])
max_mt       <- as.numeric(args[8])
max_ribo     <- as.numeric(args[9])

# Load Seurat object with doublets removed
obj <- readRDS(input)
n_cells_before <- ncol(obj)
n_features_before <- nrow(obj)

# Calculate QC metrics for filtering
obj[["percent.mt"]]     <- PercentageFeatureSet(obj, pattern = "^(COX|ND|CYTB|ATP)")
obj[["percent.ribo"]]   <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
obj[["complexity"]]     <- log10(obj$nCount_RNA / obj$nFeature_RNA)

qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "complexity")
df <- as_tibble(obj[[]], rownames="Cell.Barcode")

# Define filtering thresholds
thresholds <- list(
  nCount_RNA = c(min = min_counts, max = max_counts),
  nFeature_RNA = c(min = min_features, max = max_features),
  percent.mt = c(max = max_mt),
  percent.ribo = c(max = max_ribo)
)

# Generate QC plots with threshold lines
pdf(file.path(outdir, paste0(id, "_plots.pdf")), width = 8, height = 6)

for (dir in qc_metrics) {
  dir.create(file.path(outdir, dir), recursive = TRUE, showWarnings = FALSE)
}

# QC plots
for (i in qc_metrics) {
  # Define axis labels and units
  y_label <- switch(i,
    "nCount_RNA" = "UMI Count (counts)",
    "nFeature_RNA" = "Feature Count (genes)",
    "percent.mt" = "Mitochondrial Content (%)",
    "percent.ribo" = "Ribosomal Content (%)",
    "complexity" = "Complexity (log10)"
  )

  plot_title <- switch(i,
    "nCount_RNA" = "UMI Count per Cell Distribution",
    "nFeature_RNA" = "Feature Count per Cell Distribution",
    "percent.mt" = "Mitochondrial Percentage Distribution",
    "percent.ribo" = "Ribosomal Percentage Distribution",
    "complexity" = "Library Complexity Distribution"
  )

  # violin plots
  p <- VlnPlot(obj, features = i, layer = "counts", pt.size=0.05) +
    labs(
      x = NULL,
      y = y_label,
      title = paste(plot_title, "-", id)
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.position = "none"
    )

  # Add threshold lines
  if (i %in% names(thresholds)) {
    if ("min" %in% names(thresholds[[i]])) {
      p <- p + geom_hline(yintercept = thresholds[[i]]["min"],
                         color = "red", linetype = "dashed", size = 0.8)
    }
    if ("max" %in% names(thresholds[[i]])) {
      p <- p + geom_hline(yintercept = thresholds[[i]]["max"],
                         color = "red", linetype = "dashed", size = 0.8)
    }
  }

  ggsave(file.path(outdir, i, paste0(i, "_vln_", id, ".png")), p,
         width = 10.67, height = 8, dpi = 150, bg = "white")
  print(p)


  # histograms
  p <- ggplot(obj@meta.data, aes_string(x = i)) +
    geom_histogram(bins = 100, fill = "gray70", color = "black") +
    labs(
      x = y_label,
      y = "Number of Cells (count)",
      title = paste(plot_title, "Histogram -", id)
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14)
    )

  # Add threshold lines
  if (i %in% names(thresholds)) {
    if ("min" %in% names(thresholds[[i]])) {
      p <- p + geom_vline(xintercept = thresholds[[i]]["min"],
                         color = "red", linetype = "dashed", size = 0.8)
    }
    if ("max" %in% names(thresholds[[i]])) {
      p <- p + geom_vline(xintercept = thresholds[[i]]["max"],
                         color = "red", linetype = "dashed", size = 0.8)
    }
  }

  ggsave(file.path(outdir, i, paste0(i, "_dist_", id, ".png")), p,
         width = 10.67, height = 8, dpi = 150, bg = "white")
  print(p)


  # kernel density estimate
  p <- ggplot(obj@meta.data, aes_string(x = i)) +
    geom_density(fill = "skyblue", alpha = 0.6) +
    scale_x_log10() +
    labs(
      x = paste(y_label, "(log10 scale)"),
      y = "Density",
      title = paste(plot_title, "Kernel Density Estimate -", id)
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14)
    )

  # Add threshold lines
  if (i %in% names(thresholds)) {
    if ("min" %in% names(thresholds[[i]])) {
      p <- p + geom_vline(xintercept = thresholds[[i]]["min"],
                         color = "red", linetype = "dashed", size = 0.8)
    }
    if ("max" %in% names(thresholds[[i]])) {
      p <- p + geom_vline(xintercept = thresholds[[i]]["max"],
                         color = "red", linetype = "dashed", size = 0.8)
    }
  }

  ggsave(file.path(outdir, i, paste0(i, "_KDE_", id, ".png")), p,
         width = 10.67, height = 8, dpi = 150, bg = "white")
  print(p)


  # scatter plots
  if (i != "nCount_RNA") {
    y_label_scatter <- y_label

    # Calculate correlation
    cor_value <- cor(obj@meta.data$nCount_RNA, obj@meta.data[[i]], use = "complete.obs")

    p <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = i, pt.size = 0.5) +
      labs(
        x = "UMI Count (counts)",
        y = y_label_scatter,
        title = paste(y_label_scatter, "vs UMI Count \nr =", round(cor_value, 3), "-", id),
        color = "Cell Density"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.position = "none"
      )

    # Add threshold lines
    if (i %in% names(thresholds)) {
      if ("min" %in% names(thresholds[[i]])) {
        p <- p + geom_hline(yintercept = thresholds[[i]]["min"],
                           color = "red", linetype = "dashed", size = 0.8)
      }
      if ("max" %in% names(thresholds[[i]])) {
        p <- p + geom_hline(yintercept = thresholds[[i]]["max"],
                           color = "red", linetype = "dashed", size = 0.8)
      }
    }
    # Add nCount_RNA thresholds as vertical lines
    p <- p +
      geom_vline(xintercept = thresholds[["nCount_RNA"]]["min"],
                 color = "red", linetype = "dashed", size = 0.8) +
      geom_vline(xintercept = thresholds[["nCount_RNA"]]["max"],
                 color = "red", linetype = "dashed", size = 0.8)

    ggsave(file.path(outdir, i, paste0(i, "_vs_nCount_RNA_", id, ".png")), p,
           width = 10.67, height = 8, dpi = 150, bg = "white")
    print(p)
  }


  # scatter plots with mt and ribo
  if (i %in% c("percent.mt","percent.ribo")) {
    color_label <- switch(i,
      "percent.mt" = "Mitochondrial\nContent (%)",
      "percent.ribo" = "Ribosomal\nContent (%)"
    )

    plot_title_colored <- switch(i,
      "percent.mt" = "Feature vs UMI Count colored by Mitochondrial Content",
      "percent.ribo" = "Feature vs UMI Count colored by Ribosomal Content"
    )

    # Calculate correlation between nCount_RNA and nFeature_RNA
    cor_value <- cor(obj@meta.data$nCount_RNA, obj@meta.data$nFeature_RNA, use = "complete.obs")

    p <- df %>% arrange(!!sym(i)) %>%
      ggplot(aes_string("nCount_RNA", "nFeature_RNA", colour = i)) +
      geom_point(size = 1, alpha = 0.6) +
      scale_color_gradientn(
        colors = c("black","blue","green2","red","yellow"),
        name = color_label,
        guide = guide_colorbar(barwidth = 1.5, barheight = 15)
      ) +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        x = "UMI Count (counts, log10 scale)",
        y = "Feature Count (genes, log10 scale)",
        title = paste(plot_title_colored, "\nr =", round(cor_value, 3), "-", id)
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
      ) +
      # Add UMI count thresholds (vertical lines)
      geom_vline(xintercept = thresholds[["nCount_RNA"]]["min"],
                 color = "red", linetype = "dashed", size = 0.8) +
      geom_vline(xintercept = thresholds[["nCount_RNA"]]["max"],
                 color = "red", linetype = "dashed", size = 0.8) +
      # Add feature count thresholds (horizontal lines)
      geom_hline(yintercept = thresholds[["nFeature_RNA"]]["min"],
                 color = "red", linetype = "dashed", size = 0.8) +
      geom_hline(yintercept = thresholds[["nFeature_RNA"]]["max"],
                 color = "red", linetype = "dashed", size = 0.8)

    ggsave(file.path(outdir, i, paste0(i, "_feature_vs_count_", id, ".png")), p,
           width = 10.67, height = 8, dpi = 150, bg = "white")
    print(p)
  }
}

dev.off()

# Filter
obj <- subset(
  obj,
  subset =
    nCount_RNA   >= min_counts   &
    nCount_RNA   <= max_counts   &
    nFeature_RNA >= min_features &
    nFeature_RNA <= max_features &
    percent.mt   <= max_mt       &
    percent.ribo <= max_ribo
)

# Save object
saveRDS(obj, file = file.path(dirname(outdir), paste0(id, ".rds")))

# Summary table
write.table(
  data.frame(
    id,
    n_cells_before,
    n_cells_after = ncol(obj),
    n_features_before,
    n_features_after = nrow(obj),
    min_counts,
    max_counts,
    min_features,
    max_features,
    max_mt,
    max_ribo
  ),
  file = file.path(outdir, paste0(id, "_filtering_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
