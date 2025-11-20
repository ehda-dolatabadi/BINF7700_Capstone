#!/bin/bash
#SBATCH --job-name=cluster
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --chdir=/scratch/dolatabadi.e   # log directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log

set -euo pipefail

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Running PCA..."
Rscript "$WORK/scripts/seurat/05_pca.R"			"$ID"	"$OUT"	"$OUT/${ID}_04_integrated.rds"
echo "PCA completed✅"

echo "Clustering..."
Rscript "$WORK/scripts/seurat/06_cluster.R"           	"$ID"	"$OUT"	"$OUT/${ID}_05_pca.rds"
echo "Clustering completed✅"

echo "Running UMAP..."
Rscript "$WORK/scripts/seurat/07_umap.R"               	"$ID"	"$OUT"	"$OUT/${ID}_06_clustered.rds"
echo "UMAP completed✅"

echo "Finding Markers..."
Rscript "$WORK/scripts/seurat/08_find_markers.R"        "$ID"	"$OUT"	"$OUT/${ID}_06_clustered.rds"
echo "Markers finding completed✅"

echo "Scoring markers..."
Rscript "$WORK/scripts/seurat/09_score_markers.R"	"$ID"   "$OUT"  "$OUT/${ID}_06_clustered.rds"
echo "Scoring completed✅"

# Sync
rsync -av --update --info=progress2 --stats --exclude '*.rds' "$OUT/" "$WORK/results/seurat/"
