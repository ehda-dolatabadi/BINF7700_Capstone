#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --chdir=/scratch/dolatabadi.e   # log directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log
#SBATCH --array=0-7

set -euo pipefail

# Variables
SMP="${SMPS[$SLURM_ARRAY_TASK_ID]}"

DIR="$RESULT/cellranger/$SMP/outs/filtered_feature_bc_matrix"
OUT_map="$OUT/${ID}_00_mapped"
OUT_qc="$OUT/${ID}_01_qc/$SMP"
OUT_filter="$OUT/${ID}_02_filtered/$SMP"
OUT_normalize="$OUT/${ID}_03_normalized"

mkdir -p "$OUT_map" "$OUT_qc" "$OUT_filter" "$OUT_normalize"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Workflow
echo "Mapping features of $SMP..."
Rscript "$WORK/scripts/seurat/00_map_features.R"	"$SMP"	"$OUT_map"	"$DIR"	"$TSV"
echo "Mapping of $SMP completed✅"

echo "Running QC of $SMP..."
Rscript "$WORK/scripts/seurat/01_qc.R"                  "$SMP" "$OUT_qc"	"$OUT_map/${SMP}_mapped.rds"
echo "QC of $SMP completed✅"

echo "Filtering of $SMP..."
Rscript "$WORK/scripts/seurat/02_filter.R"              "$SMP" "$OUT_filter"	"$OUT_map/${SMP}_mapped.rds"
echo "Filtering of $SMP completed✅"

echo "Normalizing of $SMP..."
Rscript "$WORK/scripts/seurat/03_normalize.R"           "$SMP" "$OUT_normalize"	"$OUT_filter/${SMP}_filtered.rds"
echo "Normalization of $SMP completed✅"
