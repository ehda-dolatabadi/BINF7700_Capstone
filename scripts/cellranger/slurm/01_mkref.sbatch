#!/bin/bash
#SBATCH --job-name=mkref
#SBATCH --partition=courses
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --array=0-1

set -euo pipefail

# References
FA_names=("GCF_040938575.1_UKY_AmexF1_1_genomic.fna" "AmexG_v6.0-DD.fa")
GTF_names=("GCF_040938575.1_UKY_AmexF1_1_genomic.gtf" "AmexT_v47-AmexG_v6.0-DD.noWS.gtf")
REF_names=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")

# Paths
source "../../../config/default_paths.sh"
if [ -f "../../../config/local_paths.sh" ]; then
    source "../../../config/local_paths.sh"
fi

OUT_ref="$DATA/ref_genomes/cellranger_mkref"

# Derive indices
FA="$REF_FILES/${FA_names[$SLURM_ARRAY_TASK_ID]}"
GTF="$REF_FILES/${GTF_names[$SLURM_ARRAY_TASK_ID]}"
REF="${REF_names[$SLURM_ARRAY_TASK_ID]}"

mkdir -p "$OUT_ref"
cd "$OUT_ref"

rm -rf ./"$REF"/
rm -rf ./"mkref_${REF}"/

# Run pipeline
cellranger mkref \
  --fasta="$FA" \
  --genes="$GTF" \
  --genome="$REF" \
  --nthreads="$SLURM_CPUS_PER_TASK" \
  --memgb=$(( SLURM_MEM_PER_NODE / 1024 - 20))

rm -r ./"mkref_${REF}"/
