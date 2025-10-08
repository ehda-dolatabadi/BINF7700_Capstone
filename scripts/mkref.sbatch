#!/bin/bash
#SBATCH --job-name=mkref
#SBATCH --partition=courses
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --chdir=/scratch/dolatabadi.e   # working directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log

set -euo pipefail

# Variables
SRC="$WORK/ref_genomes/ref_files"

FA="$SRC/GCF_040938575.1_UKY_AmexF1_1_genomic.fna"
GTF="$SRC/GCF_040938575.1_UKY_AmexF1_1_genomic.gtf"

GNM_name="UKY_AmexF1_1_genomic_$(date +%Y%m%d)_${SLURM_JOB_ID}"

# Paths
FA_name=$(basename "$FA")
GTF_name=$(basename "$GTF")

OUT="$WORK/ref_genomes/cellranger_mkref/$GNM_name"

TMP="$SCRATCH/ref_genomes"
FILES="$TMP/ref_files"

# Sync files to scratch
echo "Syncing files to scratch..."

mkdir -p "$FILES"
cp -rn "$FA" "$GTF" "$FILES"/

cd "$TMP"
rm -rf ./"$GNM_name"/
rm -rf ./"mkref_${GNM_name}"/

# Run pipeline
# cpu defaults to 1
# memgb REQUIRED. defaults to 16. reserve ~16 GB headroom for OS
cellranger mkref \
  --fasta="$FILES/$FA_name" \
  --genes="$FILES/$GTF_name" \
  --genome="$GNM_name" \
  --nthreads="$SLURM_CPUS_PER_TASK" \
  --memgb=$(( SLURM_MEM_PER_NODE / 1024 - 16 ))

rsync -a --info=progress2 "$GNM_name"/ "$OUT"/
rm -r ./"mkref_${GNM_name}"/
