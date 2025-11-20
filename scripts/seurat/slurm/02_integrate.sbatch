#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --chdir=/scratch/dolatabadi.e   # log directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log

set -euo pipefail

# Variables
OUT_normalize="$OUT/${ID}_03_normalized"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Integrting normalized samples..."
Rscript "$WORK/scripts/seurat/04_integrate.R" "$ID" "$OUT" "$OUT_normalize"/*_normalized.rds
echo "Integration completed✅"
