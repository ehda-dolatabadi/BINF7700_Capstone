#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH --partition=courses
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --chdir=/scratch/dolatabadi.e   # working directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log

set -euo pipefail

# Variables
REF_source="$WORK/ref_genomes/cellranger_mkref/AmexT_v47-AmexG_v6_0-DD"
FQ_source="$WORK/raw/Li_dataset/7DPA"

# Sync files and change directory to scratch
REF_name=$(basename "$REF_source")
#FQ_name=$(basename "$FQ_source")
FQ_name="7dpa"

SMP="$FQ_name"
REF="$SCRATCH/ref_genomes/cellranger_mkref/$REF_name"
FQ="$SCRATCH/fastqs/$FQ_name"
TMP="$SCRATCH/outputs/cellranger/${REF_name}_${FQ_name}_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID}"
OUT="$WORK/outputs/cellranger/${REF_name}_${FQ_name}_$(date +%Y%m%d)_${SLURM_JOB_ID}"

mkdir -p "$REF" "$FQ" "$TMP" "$OUT"
cd "$TMP"

rsync -a --ignore-existing --info=progress2 "$REF_source"/ "$REF"/
rsync -a --ignore-existing --info=progress2 "$FQ_source"/ "$FQ"/

# Run pipeline
cellranger count \
  --id=test_count \
  --transcriptome="$REF" \
  --fastqs="$FQ" \
  --sample="$SMP" \
  --create-bam=false \
  --nosecondary \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem="$(( SLURM_MEM_PER_NODE / 1024 - 16 ))"		# reserve ~16 GB headroom for OS
  # --chemistry=SC3Pv3 is required if dataset is old versions

rsync -a --ignore-existing --info=progress2 "$TMP"/ "$OUT"/
