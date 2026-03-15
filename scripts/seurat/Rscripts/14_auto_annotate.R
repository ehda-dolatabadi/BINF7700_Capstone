# 1. Install dependencies ─────────────────────────────────────────────────────

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleR", "celldex", "scuttle"))

library(Seurat)
library(SingleR)
library(celldex)
library(scuttle)
library(ggplot2)


# 2. Load reference dataset ───────────────────────────────────────────────────
# Choose the reference that matches your organism and tissue type:
#   HumanPrimaryCellAtlasData()  — broad human cell types
#   BlueprintEncodeData()        — human, good for immune and stromal
#   MouseRNAseqData()            — mouse bulk RNA reference
#   ImmGenData()                 — mouse immune cells (high resolution)
#   DatabaseImmuneCellExpressionData() — human immune cells

ref <- celldex::HumanPrimaryCellAtlasData()


# 3. Prepare query matrix ─────────────────────────────────────────────────────

query <- GetAssayData(seu, assay = "RNA", slot = "data")


# 4. Run SingleR ──────────────────────────────────────────────────────────────
# label.main = broad cell types
# label.fine = finer subtypes (swap if higher resolution is needed)
# fine.tune  = TRUE improves accuracy on ambiguous cells

pred <- SingleR(
  test      = query,
  ref       = ref,
  labels    = ref$label.main,
  fine.tune = TRUE
)

table(pred$labels)


# 5. Add labels to Seurat object ──────────────────────────────────────────────

seu$singler_label        <- pred$labels
seu$singler_label_pruned <- pred$pruned.labels   # NA = ambiguous / low confidence
seu$singler_score        <- apply(pred$scores, 1, max)

cat("Total cells annotated:", sum(!is.na(pred$pruned.labels)), "\n")
cat("Ambiguous cells (NA):", sum(is.na(pred$pruned.labels)), "\n")


# 6. Visualize ────────────────────────────────────────────────────────────────

# UMAP colored by annotation
DimPlot(
  seu,
  reduction = "umap",
  group.by  = "singler_label",
  label     = TRUE,
  repel     = TRUE,
  pt.size   = 0.3
) +
  ggtitle("SingleR Cell Type Annotation") +
  theme(plot.title = element_text(hjust = 0.5))

# Score distribution per cell type
ggplot(data.frame(label = pred$labels, score = apply(pred$scores, 1, max)),
       aes(x = label, y = score)) +
  geom_violin(fill = "steelblue", alpha = 0.7) +
  geom_jitter(size = 0.1, alpha = 0.2) +
  coord_flip() +
  labs(
    title = "Annotation Confidence by Cell Type",
    x     = "Cell Type",
    y     = "Score"
  ) +
  theme_classic()

# Heatmap of scores across all cell types
plotScoreHeatmap(pred)

# Delta distribution — low delta = uncertain assignment
plotDeltaDistribution(pred, ncol = 3)
