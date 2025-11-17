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
id	<- args[1]
outdir	<- args[2]
input	<- args[3]

# Parameters
epi_thr <- 1.5
ery_thr <- 10

# Load clustered object
obj <- readRDS(input)

# Define your custom order
timepoint_order <- c("control", "3h", "24h", "72h", "7dpa", "14dpa", "22dpa", "33dpa")
obj$orig.ident <- factor(obj$orig.ident, levels = timepoint_order)

# Define marker sets
epithelial_markers <- c("KRT8", "KRT18", "KRT1", "EPCAM", "CLDN7", 
                        "CDH1", "DSP", "DSG2", "EPPK1", "S100P")

erythrocyte_markers <- c("HBA1", "HBB", "HBA2", "HEMGN", 
                         "ALAS2", "PRDX2", "ANK1", "HBBP1")

schwann_markers <- c("SOX10", "S100", "S100B", "NGFR", "p75NTR", "MPZ", "MBP", "PMP22", "PLP1", "PRX", 
		     "NCAM", "NCAM1", "L1CAM", "SCN7A", "SOX2", "GAP43", "EGR2", "Krox20", "POU3F1", "OCT6")

epithelial_markers_available <- epithelial_markers[epithelial_markers %in% rownames(obj)]
erythrocyte_markers_available <- erythrocyte_markers[erythrocyte_markers %in% rownames(obj)]
schwann_markers_available <- schwann_markers[schwann_markers %in% rownames(obj)]

# Add module score for epithelial cells
obj <- AddModuleScore(
  obj,
  features = list(epithelial_markers_available),
  name = "Epithelial_score"
)

# Add module score for erythrocyte cells
obj <- AddModuleScore(
  obj,
  features = list(erythrocyte_markers_available),
  name = "Erythrocyte_score"
)

# Add module score for schwann cells
obj <- AddModuleScore(
  obj,
  features = list(schwann_markers_available),
  name = "Schwann_score"
)

# Plot VlnPlot
png(file.path(outdir, paste0(id, "_violin_score.png")), width = 1600, height = 1200)
VlnPlot(obj, features = c("Epithelial_score1", "Erythrocyte_score1", "Schwann_score1"), pt.size = 0, group.by = "orig.ident") +
  labs(x = NULL)
dev.off()

# Visualize distribution
png(file.path(outdir, paste0(id, "_epithelial_score_distribution.png")), width = 1600, height = 1200)
hist(obj$Epithelial_score1, breaks = 50, main = "Epithelial Score Distribution", xlab = "Epithelial Score")
abline(v = epi_thr, col = "red", lwd = 2, lty = 2)
dev.off()

png(file.path(outdir, paste0(id, "_erythrocyte_score_distribution.png")), width = 1600, height = 1200)
hist(obj$Erythrocyte_score1, breaks = 50, main = "Erythrocyte Score Distribution", xlab = "Erythrocyte Score")
abline(v = 10.0, col = "red", lwd = 2, lty = 2)
dev.off()

png(file.path(outdir, paste0(id, "_schwann_score_distribution.png")), width = 1600, height = 1200)
hist(obj$Schwann_score1, breaks = 50, main = "Schwann Score Distribution", xlab = "Schwann Score")
dev.off()

df %>% arrange(obj$Schwann_score1) %>%
      ggplot(aes_string(obj$Epithelial_score1, obj$Erythrocyte_score1, colour = obj$Schwann_score1))) +
      geom_point() +
      scale_color_gradientn(colors = c("black","blue","green2","red","yellow"),
                            guide = guide_colorbar(barwidth = 2, barheight = 30)) +
      scale_x_log10() +
      scale_y_log10()

# Filtering
#cells_to_keep <- which(obj$Epithelial_score1 < epi_thr & 
                       obj$Erythrocyte_score1 < ery_thr)

# Save object
saveRDS(obj, file = file.path(outdir, paste0(id, "_scored.rds")))
