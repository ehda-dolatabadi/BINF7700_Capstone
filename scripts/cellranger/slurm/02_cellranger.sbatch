#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH --partition=courses
#SBATCH --time=24:00:00			# limit is 24hr
#SBATCH --cpus-per-task=4		# limit is 128/node
#SBATCH --mem=200G			# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e   # working directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log
#SBATCH --array=0-15

set -euo pipefail

# Arrays
SRC="$SCRATCH"

REF_names=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")

FQ_srcs=("0HPA" "3HPA" "1DPA" "3DPA" "7DPA" "14DPA" "22DPA" "33DPA")
FQ_names=("control" "3h" "24h" "72h" "7dpa" "14dpa" "22dpa" "33dpa")

# Derive indices
ref_idx=$((SLURM_ARRAY_TASK_ID / ${#FQ_srcs[@]}))
fq_idx=$((SLURM_ARRAY_TASK_ID % ${#FQ_srcs[@]}))

REF_name=${REF_names[$ref_idx]}
FQ_src="$SRC/Li_dataset/${FQ_srcs[$fq_idx]}"
FQ_name=${FQ_names[$fq_idx]}

# Paths
REF="$SRC/ref_genomes/cellranger_mkref/$REF_name"
OUT="$SRC/outputs/ref_${REF_name}/cellranger"

mkdir -p "$OUT"
cd "$OUT"

# Run Cell Ranger on scratch
cellranger count \
  --id="$FQ_name" \
  --transcriptome="$REF" \
  --fastqs="$FQ_src" \
  --sample="$FQ_name" \
  --create-bam=false \
  --nosecondary \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem=$(( SLURM_MEM_PER_NODE / 1024 - 20))
