#!/bin/bash
# Script: 08_auto_annotate.sbatch
# Purpose: Annotate clusters with cell type labels

#SBATCH --job-name=auto_annotate
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Parameter defaults
ID="${ID:-all}"
CELLS="${CELLS:-Neurons}"

# Load configuration and setup paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

# Define script and output directories
RUN="$WORK/scripts/seurat/Rscripts"
OUT_umap="$OUT/seurat/${ID}_08_umap"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate seurat_a

# Annotation workflow
# Rename cluster identities using auto annotation
echo "Annotating cells automatically..."
Rscript "$RUN/14_auto_annotate.R" \
  "$ID" \
  "$OUT_umap" \
  "$OUT/seurat/${ID}_processed.rds" \
  "$CELLS"
echo "Auto annotation completed✅"
