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
SRC="$SCRATCH"
FQ_src="$SRC/Li_dataset"
OUT="$DATA/outputs/fastqc"

mkdir -p $OUT
cd $OUT

# Environment
eval "$(conda shell.bash hook)"
conda activate preprocess

# Run FastQC
fastqc $FQ_src/*/*.fastq.gz \
  --outdir=$OUT \
  --threads=16

# Run MultiQC
multiqc $OUT\
  --outdir $OUT

# Sync
rsync -av --update --info=progress2 --stats --exclude 'cellranger/' --exclude '*.rds' "$DATA/outputs/" "$WORK/results/"
