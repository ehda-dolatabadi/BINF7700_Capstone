#!/bin/bash
#SBATCH --job-name=mkref
#SBATCH --partition=courses
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --array=0-1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}.log"
exec 2>&1

# Derive indices
FA="$REF_FILES/${FA_names[$SLURM_ARRAY_TASK_ID]}"
GTF="$REF_FILES/${GTF_names[$SLURM_ARRAY_TASK_ID]}"
REF="${REF_names[$SLURM_ARRAY_TASK_ID]}"

OUT_ref="$DATA/ref_genomes/cellranger_mkref"
mkdir -p "$OUT_ref"

rm -rf ./"$REF"/
rm -rf ./"mkref_${REF}"/

# Run pipeline
cellranger mkref \
  --output-dir="$OUT_ref/$REF" \
  --fasta="$FA" \
  --genes="$GTF" \
  --genome="$REF" \
  --nthreads="$SLURM_CPUS_PER_TASK" \
  --memgb=$(( SLURM_MEM_PER_NODE / 1024 - 20))

rm -r ./"mkref_${REF}"/
