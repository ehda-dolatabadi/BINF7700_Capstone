#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-7

set -euo pipefail

# Default Parameters
ID="${ID:-main}"
MIN_CELLS="${MIN_CELLS:-5}"
MIN_FEATURES="${MIN_FEATURES:-500}"
MAX_FEATURES="${MAX_FEATURES:-5000}"
MIN_COUNTS="${MIN_COUNTS:-1000}"
MAX_COUNTS="${MAX_COUNTS:-30000}"
MAX_MT="${MAX_MT:-10}"
MAX_RIBO="${MAX_RIBO:-30}"

# Paths
SMPS=("control" "3h" "24h" "72h" "7dpa" "14dpa" "22dpa" "33dpa")
SMP="${SMPS[$SLURM_ARRAY_TASK_ID]}"

DIR="$DATA/outputs/cellranger/$SMP/outs/filtered_feature_bc_matrix"
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
Rscript "$RUN/00_map_features.R" \
  "$SMP" \
  "$OUT_map" \
  "$DIR" \
  "$TSV" \
  "$MIN_CELLS" \
  "$MIN_FEATURES"
echo "Mapping of $SMP completed✅"

echo "Running QC of $SMP..."
Rscript "$RUN/01_qc.R" \
  "$SMP" \
  "$OUT_qc" \
  "$OUT_map/${SMP}_mapped.rds"
echo "QC of $SMP completed✅"

echo "Filtering of $SMP..."
Rscript "$RUN/02_filter.R" \
  "$SMP" \
  "$OUT_filter" \
  "$OUT_map/${SMP}_mapped.rds" \
  "$MIN_COUNTS" \
  "$MAX_COUNTS" \
  "$MIN_FEATURES" \
  "$MAX_FEATURES" \
  "$MAX_MT" \
  "$MAX_RIBO"
echo "Filtering of $SMP completed✅"

echo "Normalizing of $SMP..."
Rscript "$RUN/03_normalize.R" \
  "$SMP" \
  "$OUT_normalize" \
  "$OUT_filter/${SMP}_filtered.rds"
echo "Normalization of $SMP completed✅"

# Sync
rsync -av --info=progress2 --stats --exclude '*.rds' "$OUT/" "$WORK/results/seurat/"
