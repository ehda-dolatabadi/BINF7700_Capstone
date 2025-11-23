#!/bin/bash
#SBATCH --job-name=integrate
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Parameters
ID="${ID:-main}"

# Paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}.log"
exec 2>&1

RUN="$WORK/scripts/seurat/Rscripts"
OUT_normalize="$OUT/seurat/${ID}_03_normalized"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Integrating normalized samples..."
Rscript "$RUN/04_integrate.R" "$ID" "$OUT/seurat" "$OUT_normalize"/*_normalized.rds
echo "Integration completed✅"
