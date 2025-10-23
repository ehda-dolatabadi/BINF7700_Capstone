#!/bin/bash
#SBATCH --job-name=mkref
#SBATCH --partition=courses
#SBATCH --time=24:00:00			# limit is 24hr
#SBATCH --cpus-per-task=16		# limit is 58
#SBATCH --mem=60G			# limit is 185G
#SBATCH --chdir=/scratch/dolatabadi.e   # working directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log
#SBATCH --array=0-1

set -euo pipefail

# Arrays
SRC="$SCRATCH/ref_genomes/ref_files"

FA_names=("GCF_040938575.1_UKY_AmexF1_1_genomic.fna" "AmexG_v6.0-DD.fa")
GTF_names=("GCF_040938575.1_UKY_AmexF1_1_genomic.gtf" "AmexT_v47-AmexG_v6.0-DD.noWS.gtf")
GNM_names=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")

# Derive indices
FA="$SRC/${FA_names[$SLURM_ARRAY_TASK_ID]}"
GTF="$SRC/${GTF_names[$SLURM_ARRAY_TASK_ID]}"
GNM_name="${GNM_names[$SLURM_ARRAY_TASK_ID]}"

# Paths
OUT="$SCRATCH/ref_genomes/cellranger_mkref"

mkdir -p "$OUT"
cd "$OUT"

rm -rf ./"$GNM_name"/
rm -rf ./"mkref_${GNM_name}"/

# Run pipeline
cellranger mkref \
  --fasta="$FA" \
  --genes="$GTF" \
  --genome="$GNM_name" \
  --nthreads="$SLURM_CPUS_PER_TASK" \
  --memgb=$(( SLURM_MEM_PER_NODE / 1024 - 20))

rm -r ./"mkref_${GNM_name}"/
