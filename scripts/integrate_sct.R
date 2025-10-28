#!/usr/bin/env Rscript
# Purpose: Integrate multiple SCT-normalized Seurat objects

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(777)

# Args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript integrate_sct.R <project_id> <outdir> <sample1_sct.rds> <sample2_sct.rds> ...")
}

project_id <- args[1]
outdir     <- args[2]
rds_paths  <- args[3:length(args)]

# Load SCT objects
obj_list <- lapply(rds_paths, readRDS)

# Integration features
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)

# Find anchors and integrate (batch correction)
anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  normalization.method = "SCT",
  anchor.features = features
)

integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

# Use the integrated assay for downstream steps
DefaultAssay(integrated) <- "integrated"

# Save
saveRDS(integrated, file = file.path(outdir, paste0(project_id, "_integrated.rds")))

# Summary TSV
n_cells_total <- sum(sapply(obj_list, ncol))
write.table(
  data.frame(
    project_id            = project_id,
    n_samples             = length(obj_list),
    n_cells_total         = n_cells_total,
    n_features_integrated = nrow(integrated),
    default_assay         = DefaultAssay(integrated),
    n_integration_features= length(features)
  ),
  file = file.path(outdir, paste0(project_id, "_integration_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(sessionInfo())
