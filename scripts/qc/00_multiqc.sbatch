#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --partition=courses
#SBATCH --time=24:00:00			# limit is 24hr
#SBATCH --cpus-per-task=18		# limit is 128/node
#SBATCH --mem=40G			# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e   # working directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log

set -euo pipefail

# Paths
OUT_qc="$OUT/fastqc"

mkdir -p "$OUT_qc"
cd "$OUT_qc"

# Environment
eval "$(conda shell.bash hook)"
conda activate preprocess

# Run FastQC
fastqc $FQ_src/*/*.fastq.gz \
  --outdir="$OUT_qc" \
  --threads=16

# Run MultiQC
multiqc "$OUT_qc"\
  --outdir "$OUT_qc"

# Sync
rsync -av --info=progress2 --stats --exclude 'multiqc_data/' --exclude '*.zip' "$OUT_qc/" "$WORK/results/qc/"
