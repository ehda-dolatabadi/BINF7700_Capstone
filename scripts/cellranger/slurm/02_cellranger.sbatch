#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH --partition=courses
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=240G
#SBATCH --array=0-15
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Paths
cd "$SLURM_SUBMIT_DIR"
source "config/default_paths.sh"
[ -f "config/local_paths.sh" ] && source "config/local_paths.sh"

exec > "$LOG/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec 2>&1

OUT_ref="$DATA/ref_genomes/cellranger_mkref"

# Derive indices
ref_idx=$((SLURM_ARRAY_TASK_ID / ${#FQ_names[@]}))
fq_idx=$((SLURM_ARRAY_TASK_ID % ${#FQ_names[@]}))

REF_name=${REF_names[$ref_idx]}
FQ_name="$FQ_src/${FQ_names[$fq_idx]}"
SMP=${SMPS[$fq_idx]}

REF="$OUT_ref/$REF_name"
OUT_cellranger="$OUT/cellranger/ref_${REF_name}"
mkdir -p "$OUT_cellranger"
cd "$OUT_cellranger"

# Run Pipeline
cellranger count \
  --output-dir="$OUT_cellranger" \
  --id="$SMP" \
  --transcriptome="$REF" \
  --fastqs="$FQ_name" \
  --sample="$SMP" \
  --create-bam=false \
  --nosecondary \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem=$(( SLURM_MEM_PER_NODE / 1024 - 20))

rsync -av --info=progress2 --stats "$OUT_cellranger/$SMP/outs/web_summary.html" "$WORK/results/cellranger/ref_${REF_name}/$SMP/"
