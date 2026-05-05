#!/bin/bash
# Script: seurat_wrapper.sbatch
# Purpose: Install SeuratWrappers in R environment if not installed

#SBATCH --job-name=seurat_wrappers
#SBATCH --partition=courses
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate monocle3

# Install SeuratWrappers
Rscript -e '
library(remotes)
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
    remotes::install_github("satijalab/seurat-wrappers", upgrade = "never")
}
'
