#!/usr/bin/env Rscript
# Script: 01_qc.R
# Purpose: Run quality control analysis on a mapped sample
# Description: Calculates QC metrics (mitochondrial %, ribosomal %, complexity) and
#              generates comprehensive visualization plots for quality assessment
# Usage: Rscript 01_qc.R <id> <output> <input>

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
id	<- args[1]
outdir	<- args[2]
input	<- args[3]

# Load Seurat object with mapped features
obj <- readRDS(input)

# Calculate quality control metrics
obj[["percent.mt"]]	<- PercentageFeatureSet(obj, pattern = "^(COX|ND|CYTB|ATP)")
obj[["percent.ribo"]]	<- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
obj[["percent.rrna"]]	<- PercentageFeatureSet(obj, pattern = "^RRN")
obj[["percent.trna"]]	<- PercentageFeatureSet(obj, pattern = "^TRNA")
obj[["complexity"]]	<- log10(obj$nCount_RNA / obj$nFeature_RNA)

# Define QC metrics to plot
qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "complexity", "percent.rrna", "percent.trna")
df <- as_tibble(obj[[]], rownames="Cell.Barcode")

# Generate QC plots
pdf(file.path(outdir, paste0(id, "_plots.pdf")), width = 8, height = 6)

for (dir in qc_metrics) {
  dir.create(file.path(outdir, dir), recursive = TRUE, showWarnings = FALSE)
}

# Library prep control
for (i in c("percent.rrna","percent.trna")) {
  p <- VlnPlot(obj, features = i, layer = "counts", pt.size=0.05) +
    labs(
      x = NULL,
      y = "Percentage (%)",
      title = paste(gsub("\\.", " ", gsub("percent\\.", "", i)), "Distribution -", id)
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.position = "none"
    )

  ggsave(file.path(outdir, i, paste0(i, "_vln_", id, ".png")), p,
         width = 10.67, height = 8, dpi = 150, bg = "white")
  print(p)
}


# QC plots
for (i in qc_metrics[1:5]) {
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

  ggsave(file.path(outdir, i, paste0(i, "_dist_", id, ".png")), p,
         width = 10.67, height = 8, dpi = 150, bg = "white")
  print(p)
  
  
  # histogram for high counts
  if (i == "nCount_RNA") {
    p <- ggplot(subset(obj, nCount_RNA > 10000)@meta.data, aes_string(x = i)) +
      geom_histogram(bins = 100, fill = "gray70", color = "black") +
      labs(
        x = "UMI Count (counts)",
        y = "Number of Cells (count)",
        title = paste("UMI Count Distribution - High Counts (>10,000) -", id)
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)
      )

    ggsave(file.path(outdir, i, paste0(i, "_dist_high_", id, ".png")), p,
           width = 10.67, height = 8, dpi = 150, bg = "white")
    print(p)


    # histogram for low counts
    p <- ggplot(subset(obj, nCount_RNA < 2500)@meta.data, aes_string(x = i)) +
      geom_histogram(bins = 100, fill = "gray70", color = "black") +
      labs(
        x = "UMI Count (counts)",
        y = "Number of Cells (count)",
        title = paste("UMI Count Distribution - Low Counts (<2,500) -", id)
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)
      )

    ggsave(file.path(outdir, i, paste0(i, "_dist_low_", id, ".png")), p,
           width = 10.67, height = 8, dpi = 150, bg = "white")
    print(p)
  }
  
  
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

  ggsave(file.path(outdir, i, paste0(i, "_KDE_", id, ".png")), p,
         width = 10.67, height = 8, dpi = 150, bg = "white")
  print(p)
  

  # KDE for high counts
  if(i == "nCount_RNA") {
    p <- ggplot(subset(obj, nCount_RNA > 10000)@meta.data, aes_string(x = i)) +
      geom_density(fill = "skyblue", alpha = 0.6) +
      scale_x_log10() +
      labs(
        x = "UMI Count (counts, log10 scale)",
        y = "Density",
        title = paste("UMI Count KDE - High Counts (>10,000) -", id)
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)
      )

    ggsave(file.path(outdir, i, paste0(i, "_KDE_high_", id, ".png")), p,
           width = 10.67, height = 8, dpi = 150, bg = "white")
    print(p)


    # KDE for low counts
    p <- ggplot(subset(obj, nCount_RNA < 2500)@meta.data, aes_string(x = i)) +
      geom_density(fill = "skyblue", alpha = 0.6) +
      scale_x_log10() +
      labs(
        x = "UMI Count (counts, log10 scale)",
        y = "Density",
        title = paste("UMI Count KDE - Low Counts (<2,500) -", id)
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)
      )

    ggsave(file.path(outdir, i, paste0(i, "_KDE_low_", id, ".png")), p,
           width = 10.67, height = 8, dpi = 150, bg = "white")
    print(p)
  }
  
  
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
      )

    ggsave(file.path(outdir, i, paste0(i, "_feature_vs_count_", id, ".png")), p,
           width = 10.67, height = 8, dpi = 150, bg = "white")
    print(p)
  }
}

dev.off()
