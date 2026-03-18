#!/bin/bash
# Script: 01_preprocess.sbatch
# Purpose: Preprocess individual samples through the initial pipeline stages

#SBATCH --job-name=preprocess
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-7
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Parameter defaults
ID="${ID:-all}"
REF="${REF:-UKY_AmexF1_1_genomic}"
MIN_CELLS="${MIN_CELLS:-5}"
MIN_FEATURES="${MIN_FEATURES:-500}"
MAX_FEATURES="${MAX_FEATURES:-5000}"
MIN_COUNTS="${MIN_COUNTS:-1000}"
MAX_COUNTS="${MAX_COUNTS:-30000}"
MAX_MT="${MAX_MT:-10}"
MAX_RIBO="${MAX_RIBO:-35}"

# Load configuration and setup paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

# Redirect output to log file
exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec 2>&1

# Define script and reference paths
RUN="$WORK/scripts/seurat/Rscripts"
TSV="$REF_FILES/loc_map.tsv"

# Select sample based on array task ID
SMPS=("control" "3h" "24h" "72h" "7dpa" "14dpa" "22dpa" "33dpa")
SMP="${SMPS[$SLURM_ARRAY_TASK_ID]}"

# Define input and output directories
DIR="$OUT/cellranger/ref_${REF}/$SMP/outs/filtered_feature_bc_matrix"
OUT_map="$OUT/seurat/${ID}_00_mapped"
OUT_qc="$OUT/seurat/${ID}_01_qc"
OUT_db="$OUT/seurat/${ID}_02_DB_removed"
OUT_filter="$OUT/seurat/${ID}_03_filtered"
OUT_normalize="$OUT/seurat/${ID}_04_normalized"

# Create output directories
mkdir -p "$OUT_map" "$OUT_qc" "$OUT_db" "$OUT_filter" "$OUT_normalize"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Preprocessing workflow
# Step 1: Map gene IDs to gene symbols
echo "Mapping features of $SMP..."
Rscript "$RUN/00_map_features.R" \
  "$SMP" \
  "$OUT_map" \
  "$DIR" \
  "$TSV" \
  "$MIN_CELLS" \
  "$MIN_FEATURES"
echo "Mapping of $SMP completed✅"

# Step 2: Run quality control analysis
echo "Running QC of $SMP..."
Rscript "$RUN/01_qc.R" \
  "$SMP" \
  "$OUT_qc" \
  "$OUT/seurat/${SMP}.rds"
echo "QC of $SMP completed✅"

# Step 3: Detect and remove doublets using scDblFinder
echo "Removing doublets from $SMP..."
Rscript "$RUN/02_remove_doublets.R" \
  "$SMP" \
  "$OUT_db" \
  "$OUT/seurat/${SMP}.rds"
echo "Doublet removal of $SMP completed✅"

# Step 4: Filter cells based on QC thresholds
echo "Filtering of $SMP..."
Rscript "$RUN/03_filter.R" \
  "$SMP" \
  "$OUT_filter" \
  "$OUT/seurat/${SMP}.rds" \
  "$MIN_COUNTS" \
  "$MAX_COUNTS" \
  "$MIN_FEATURES" \
  "$MAX_FEATURES" \
  "$MAX_MT" \
  "$MAX_RIBO"
echo "Filtering of $SMP completed✅"

# Step 5: Normalize data using SCTransform
echo "Normalizing of $SMP..."
Rscript "$RUN/04_normalize.R" \
  "$SMP" \
  "$OUT_normalize" \
  "$OUT/seurat/${SMP}.rds"
echo "Normalization of $SMP completed✅"
