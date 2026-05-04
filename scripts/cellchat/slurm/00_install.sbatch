#!/bin/bash
# Script: 00_install.sbatch
# Purpose: Install CellChat and GitHub-only dependencies into the conda environment

#SBATCH --job-name=cellchat_install
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Load configuration and setup paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log"
exec 2>&1

# Define script directory
RUN="$WORK/scripts/cellchat/Rscripts"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate cellchat

# Install CellChat and GitHub-only dependencies
echo "Installing CellChat and dependencies..."
Rscript "$RUN/install_packages.R"
echo "Installation completed✅"

