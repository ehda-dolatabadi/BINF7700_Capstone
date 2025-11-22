#!/bin/bash
#SBATCH --job-name=integrate
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

set -euo pipefail

# Parameters
ID="${ID:-main}"

# Paths
OUT_normalize="$OUT/${ID}_03_normalized"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Integrating normalized samples..."
Rscript "$RUN/04_integrate.R" "$ID" "$OUT" "$OUT_normalize"/*_normalized.rds
echo "Integration completed✅"
