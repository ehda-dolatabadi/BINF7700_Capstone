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

# Parameters
min_cells	<- 5
min_features	<- 500

# Load and create Seurat object
counts <- Read10X(input)
obj <- CreateSeuratObject(
  counts = counts,
  min.cells = min_cells,	# keep genes detected in ≥5 cells
  min.features = min_features,	# keep cells with ≥500 genes
  project = id
)

# QC metrics
obj[["percent.mt"]]	<- PercentageFeatureSet(obj, pattern = "^(COX|ND|CYTB|ATP)")
obj[["percent.ribo"]]	<- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
obj[["percent.rrna"]]	<- PercentageFeatureSet(obj, pattern = "^RRN")
obj[["percent.trna"]]	<- PercentageFeatureSet(obj, pattern = "^TRNA")
obj[["log10UMIsPerGene"]] <- log10(obj$nCount_RNA / obj$nFeature_RNA)

qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt","percent.ribo","percent.rrna","percent.trna", "UMIsPerGene")
df <- as_tibble(obj[[]], rownames="Cell.Barcode")

# Plots
pdf(file.path(outdir, paste0(id, "_plots.pdf")), width = 15, height = 9)

# Library prep control
for (i in c("percent.rrna","percent.trna")) {
  p <- VlnPlot(obj, features = i, layer = "counts") +
    labs(x = NULL)
  
  ggsave(file.path(outdir, paste0(id, "_", i, "_vln.png")), p,
         width = 15, height = 9, dpi = 300, bg = "white")
  print(p)
}


# QC plots
# complexity
p <- VlnPlot(obj, features = "log10UMIsPerGene", layer = "counts") +
  labs(x = NULL)

ggsave(file.path(outdir, paste0(id, "_log10UMIsPerGene.png")), p,
       width = 15, height = 9, dpi = 300, bg = "white")
print(p)


for (i in qc_metrics[1:4]) {
  # violin plots
  p <- VlnPlot(obj, features = i, layer = "counts") +
    labs(x = NULL)
  
  ggsave(file.path(outdir, paste0(id, "_", i, "_vln.png")), p,
         width = 15, height = 9, dpi = 300, bg = "white")
  print(p)
  
  
  # histograms
  p <- ggplot(obj@meta.data, aes_string(x = i)) +
    geom_histogram(bins = 100, fill = "gray70", color = "black")
  
  ggsave(file.path(outdir, paste0(id, "_", i, "_dist.png")), p,
         width = 15, height = 9, dpi = 300, bg = "white")
  print(p)
  
  
  # histogram for high counts
  if (i == "nCount_RNA") {
    p <- ggplot(subset(obj, nCount_RNA > 10000)@meta.data, aes_string(x = i)) +
      geom_histogram(bins = 100, fill = "gray70", color = "black")
    
    ggsave(file.path(outdir, paste0(id, "_", i, "_dist_high.png")), p,
           width = 15, height = 9, dpi = 300, bg = "white")
    print(p)
    
    
    # histogram for low counts
    p <- ggplot(subset(obj, nCount_RNA < 2500)@meta.data, aes_string(x = i)) +
      geom_histogram(bins = 100, fill = "gray70", color = "black")
    
    ggsave(file.path(outdir, paste0(id, "_", i, "_dist_low.png")), p,
           width = 15, height = 9, dpi = 300, bg = "white")
    print(p)
  }
  
  
  # kernel density estimate
  p <- ggplot(obj@meta.data, aes_string(x = i)) +
    geom_density(fill = "skyblue") +
    scale_x_log10()
  
  ggsave(file.path(outdir, paste0(id, "_", i, "_KDE.png")), p,
         width = 15, height = 9, dpi = 300, bg = "white")
  print(p)
  

  # KDE for high counts
  if(i == "nCount_RNA") {
    p <- ggplot(subset(obj, nCount_RNA > 10000)@meta.data, aes_string(x = i)) +
      geom_density(fill = "skyblue") +
      scale_x_log10()
    
    ggsave(file.path(outdir, paste0(id, "_", i, "_KDE_high.png")), p,
           width = 15, height = 9, dpi = 300, bg = "white")
    print(p)
    
    
    # KDE for low counts
    p <- ggplot(subset(obj, nCount_RNA < 2500)@meta.data, aes_string(x = i)) +
      geom_density(fill = "skyblue") +
      scale_x_log10()
    
    ggsave(file.path(outdir, paste0(id, "_", i, "_KDE_low.png")), p,
           width = 15, height = 9, dpi = 300, bg = "white")
    print(p)
  }
  
  
  # scatter plots
  if (i != "nCount_RNA") {
    p <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = i) +
      theme(legend.title = element_blank())
    
    ggsave(file.path(outdir, paste0(id, "_", i, "_vs_nCount_RNA.png")), p,
           width = 15, height = 9, dpi = 300, bg = "white")
    print(p)
  }
  
  
  # scatter plots with mt and ribo
  if (i %in% c("percent.mt","percent.ribo")) {
    p <- df %>% arrange(i) %>% 
      ggplot(aes_string("nCount_RNA", "nFeature_RNA", colour = i)) +
      geom_point() +
      scale_color_gradientn(colors = c("black","blue","green2","red","yellow"),
                            guide = guide_colorbar(barwidth = 2, barheight = 30)) +
      scale_x_log10() + 
      scale_y_log10()
    
    ggsave(file.path(outdir, paste0(id, "_", i, "_percent.mt_feature_vs_count.png")), p,
           width = 15, height = 9, dpi = 300, bg = "white")
    print(p)
  }
}

dev.off()
