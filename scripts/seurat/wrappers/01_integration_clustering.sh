#!/bin/bash

set -euo pipefail

echo "Integrting normalized samples..." &&
#Rscript "$WORK/scripts/seurat/04_integrate.R"		"04" "$OUT" "$OUT_normalize"/*_normalized.rds &&
echo "Integration completed✅" &&

echo "Running PCA..." &&
#Rscript "$WORK/scripts/seurat/05_pca.R"			"05" "$OUT" "$OUT/04_integrated.rds" &&
echo "PCA completed✅" &&

echo "Clustering..." &&
#Rscript "$WORK/scripts/seurat/06_cluster.R"		"06" "$OUT" "$OUT/05_pca.rds" &&
echo "Clustering completed✅" &&

echo "Running UMAP..." &&
#Rscript "$WORK/scripts/seurat/07_umap.R"		"07" "$OUT" "$OUT/06_clustered.rds" &&
echo "UMAP completed✅" &&

echo "Finding Markers..." &&
Rscript "$WORK/scripts/seurat/08_find_markers.R" 	"08" "$OUT" "$OUT/06_clustered.rds" &&
echo "Markers finding completed✅"
