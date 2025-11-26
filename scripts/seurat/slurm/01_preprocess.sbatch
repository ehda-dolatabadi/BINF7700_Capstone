#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-7
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Default Parameters
ID="${ID:-all}"
REF="${REF:-UKY_AmexF1_1_genomic}"
MIN_CELLS="${MIN_CELLS:-5}"
MIN_FEATURES="${MIN_FEATURES:-500}"
MAX_FEATURES="${MAX_FEATURES:-5000}"
MIN_COUNTS="${MIN_COUNTS:-1000}"
MAX_COUNTS="${MAX_COUNTS:-30000}"
MAX_MT="${MAX_MT:-10}"
MAX_RIBO="${MAX_RIBO:-30}"

# Paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec 2>&1

RUN="$WORK/scripts/seurat/Rscripts"
TSV="$REF_FILES/loc_map.tsv"

SMPS=("control" "3h" "24h" "72h" "7dpa" "14dpa" "22dpa" "33dpa")
SMP="${SMPS[$SLURM_ARRAY_TASK_ID]}"

DIR="$OUT/cellranger/ref_${REF}/$SMP/outs/filtered_feature_bc_matrix"
OUT_map="$OUT/seurat/${ID}_00_mapped"
OUT_qc="$OUT/seurat/${ID}_01_qc/$SMP"
OUT_db="$OUT/seurat/${ID}_02_DB_removed/$SMP"
OUT_filter="$OUT/seurat/${ID}_03_filtered/$SMP"
OUT_normalize="$OUT/seurat/${ID}_04_normalized"

mkdir -p "$OUT_map" "$OUT_qc" "$OUT_db" "$OUT_filter" "$OUT_normalize"

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

echo "Removing doublets from $SMP..."
Rscript "$RUN/02_remove_doublets.R" \
  "$SMP" \
  "$OUT_db" \
  "$OUT_map/${SMP}_mapped.rds"
echo "Doublet removal of $SMP completed✅"

echo "Filtering of $SMP..."
Rscript "$RUN/03_filter.R" \
  "$SMP" \
  "$OUT_filter" \
  "$OUT_db/${SMP}_DB_removed.rds" \
  "$MIN_COUNTS" \
  "$MAX_COUNTS" \
  "$MIN_FEATURES" \
  "$MAX_FEATURES" \
  "$MAX_MT" \
  "$MAX_RIBO"
echo "Filtering of $SMP completed✅"

echo "Normalizing of $SMP..."
Rscript "$RUN/04_normalize.R" \
  "$SMP" \
  "$OUT_normalize" \
  "$OUT_filter/${SMP}_filtered.rds"
echo "Normalization of $SMP completed✅"

# Create PDF
echo "Creating PDF for $SMP..."
metrics=("nFeature_RNA" "percent.mt" "percent.ribo")
files=()
for metric in "${metrics[@]}"; do
  files+=("$OUT_qc/${SMP}_${metric}_vs_nCount_RNA.png")
  files+=("$OUT_db/${SMP}_${metric}_vs_nCount_RNA_DB_removed.png")
  files+=("$OUT_filter/${SMP}_${metric}_vs_nCount_RNA_filtered.png")
done
convert "${files[@]}" "$OUT_filter/$SMP/${SMP}_preprocessing_summary.pdf"
echo "PDF for $SMP created✅"

# Sync
rsync -av --info=progress2 --stats --exclude '*.rds' "$OUT/seurat/" "$WORK/results/seurat/"
