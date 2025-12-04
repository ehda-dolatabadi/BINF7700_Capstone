#!/usr/bin/env Rscript
# Purpose: Run QC on a 10x sample
# Usage: Rscript 1_qc.R <id> <output> <input>

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
id	<- args[1]
outdir	<- args[2]
input	<- args[3]

# Load mapped object
obj <- readRDS(input)

# QC metrics
obj[["percent.mt"]]	<- PercentageFeatureSet(obj, pattern = "^(COX|ND|CYTB|ATP)")
obj[["percent.ribo"]]	<- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
obj[["percent.rrna"]]	<- PercentageFeatureSet(obj, pattern = "^RRN")
obj[["percent.trna"]]	<- PercentageFeatureSet(obj, pattern = "^TRNA")
obj[["complexity"]]	<- log10(obj$nCount_RNA / obj$nFeature_RNA)

qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "complexity", "percent.rrna", "percent.trna")
df <- as_tibble(obj[[]], rownames="Cell.Barcode")

# Plots
pdf(file.path(outdir, paste0(id, "_plots.pdf")), width = 8, height = 6)

for (dir in qc_metrics) {
  dir.create(file.path(outdir, dir), recursive = TRUE, showWarnings = FALSE)
}

# Define thresholds
thresholds <- list(
  nCount_RNA = c(min = 1000, max = 30000),
  nFeature_RNA = c(min = 500, max = 5000),
  percent.mt = c(max = 10),
  percent.ribo = c(max = 35)
)

# Library prep control
for (i in c("percent.rrna","percent.trna")) {
  p <- VlnPlot(obj, features = i, layer = "counts") +
    labs(
      x = "Sample",
      y = "Percentage (%)",
      title = paste(gsub("\\.", " ", gsub("percent\\.", "", i)), "Distribution")
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )

  ggsave(file.path(outdir, i, paste0(i, "_vln_", id, ".png")), p,
         width = 8, height = 6, dpi = 150, bg = "white")
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
  p <- VlnPlot(obj, features = i, layer = "counts") +
    labs(
      x = "Sample",
      y = y_label,
      title = plot_title
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
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
         width = 8, height = 6, dpi = 150, bg = "white")
  print(p)
  
  
  # histograms
  p <- ggplot(obj@meta.data, aes_string(x = i)) +
    geom_histogram(bins = 100, fill = "gray70", color = "black") +
    labs(
      x = y_label,
      y = "Number of Cells (count)",
      title = paste(plot_title, "- Histogram")
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
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
         width = 8, height = 6, dpi = 150, bg = "white")
  print(p)
  
  
  # histogram for high counts
  if (i == "nCount_RNA") {
    p <- ggplot(subset(obj, nCount_RNA > 10000)@meta.data, aes_string(x = i)) +
      geom_histogram(bins = 100, fill = "gray70", color = "black") +
      labs(
        x = "UMI Count (counts)",
        y = "Number of Cells (count)",
        title = "UMI Count Distribution - High Counts (>10,000)"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      geom_vline(xintercept = thresholds[["nCount_RNA"]]["max"],
                 color = "red", linetype = "dashed", size = 0.8)

    ggsave(file.path(outdir, i, paste0(i, "_dist_high_", id, ".png")), p,
           width = 8, height = 6, dpi = 150, bg = "white")
    print(p)


    # histogram for low counts
    p <- ggplot(subset(obj, nCount_RNA < 2500)@meta.data, aes_string(x = i)) +
      geom_histogram(bins = 100, fill = "gray70", color = "black") +
      labs(
        x = "UMI Count (counts)",
        y = "Number of Cells (count)",
        title = "UMI Count Distribution - Low Counts (<2,500)"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      geom_vline(xintercept = thresholds[["nCount_RNA"]]["min"],
                 color = "red", linetype = "dashed", size = 0.8)

    ggsave(file.path(outdir, i, paste0(i, "_dist_low_", id, ".png")), p,
           width = 8, height = 6, dpi = 150, bg = "white")
    print(p)
  }
  
  
  # kernel density estimate
  p <- ggplot(obj@meta.data, aes_string(x = i)) +
    geom_density(fill = "skyblue", alpha = 0.6) +
    scale_x_log10() +
    labs(
      x = paste(y_label, "(log10 scale)"),
      y = "Density",
      title = paste(plot_title, "- Kernel Density Estimate")
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
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
         width = 8, height = 6, dpi = 150, bg = "white")
  print(p)
  

  # KDE for high counts
  if(i == "nCount_RNA") {
    p <- ggplot(subset(obj, nCount_RNA > 10000)@meta.data, aes_string(x = i)) +
      geom_density(fill = "skyblue", alpha = 0.6) +
      scale_x_log10() +
      labs(
        x = "UMI Count (counts, log10 scale)",
        y = "Density",
        title = "UMI Count KDE - High Counts (>10,000)"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      geom_vline(xintercept = thresholds[["nCount_RNA"]]["max"],
                 color = "red", linetype = "dashed", size = 0.8)

    ggsave(file.path(outdir, i, paste0(i, "_KDE_high_", id, ".png")), p,
           width = 8, height = 6, dpi = 150, bg = "white")
    print(p)


    # KDE for low counts
    p <- ggplot(subset(obj, nCount_RNA < 2500)@meta.data, aes_string(x = i)) +
      geom_density(fill = "skyblue", alpha = 0.6) +
      scale_x_log10() +
      labs(
        x = "UMI Count (counts, log10 scale)",
        y = "Density",
        title = "UMI Count KDE - Low Counts (<2,500)"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      geom_vline(xintercept = thresholds[["nCount_RNA"]]["min"],
                 color = "red", linetype = "dashed", size = 0.8)

    ggsave(file.path(outdir, i, paste0(i, "_KDE_low_", id, ".png")), p,
           width = 8, height = 6, dpi = 150, bg = "white")
    print(p)
  }
  
  
  # scatter plots
  if (i != "nCount_RNA") {
    y_label_scatter <- y_label

    p <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = i, pt.size = 0.5) +
      labs(
        x = "UMI Count (counts)",
        y = y_label_scatter,
        title = paste(y_label_scatter, "vs UMI Count"),
        color = "Cell Density"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9)
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
           width = 8, height = 6, dpi = 150, bg = "white")
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
        title = plot_title_colored
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9)
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
           width = 8, height = 6, dpi = 150, bg = "white")
    print(p)
  }
}

dev.off()
