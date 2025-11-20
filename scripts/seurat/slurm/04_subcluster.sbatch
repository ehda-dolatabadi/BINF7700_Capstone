#!/bin/bash
#SBATCH --job-name=subcluster
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

set -euo pipefail

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Subsetting data..."
Rscript "$WORK/scripts/seurat/10_subset_clusters.R"	"$ID"   "$OUT"  "$OUT/${main_ID}_06_clustered.rds" "$IDENT"
echo "Subset completed✅"

echo "Re-running normalization and integration..."
Rscript "$WORK/scripts/seurat/04_reintegrate_subset.R"  "$ID"   "$OUT"  "$OUT/${ID}_00_subset.rds"
echo "Normalization and integration completed✅"
