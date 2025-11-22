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
source "../../../config/default_paths.sh"
if [ -f "../../../config/local_paths.sh" ]; then
    source "../../../config/local_paths.sh"
fi

RUN="$WORK/scripts/seurat/Rscripts"
OUT_normalize="$OUT/seurat/${ID}_03_normalized"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Integrating normalized samples..."
Rscript "$RUN/04_integrate.R" "$ID" "$OUT/seurat" "$OUT_normalize"/*_normalized.rds
echo "Integration completed✅"
