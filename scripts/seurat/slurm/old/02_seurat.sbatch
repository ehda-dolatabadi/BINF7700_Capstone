#!/bin/bash
#SBATCH --job-name=seurat
#SBATCH --partition=courses
#SBATCH --time=24:00:00			# limit is 24hr
#SBATCH --cpus-per-task=64		# limit is 128/node
#SBATCH --mem=128G			# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e	# working directory
#SBATCH --output=logs/%x_%j.out		# stdout log
#SBATCH --error=logs/%x_%j.err		# stderr log

set -euo pipefail

# Samples
SMPS=("control" "3h" "24h" "72h" "7dpa" "14dpa" "22dpa" "33dpa")

set -a  # export all variables from here on
REF="UKY_AmexF1_1_genomic"
TSV="$DATA/ref_genomes/ref_files/loc_map.tsv"
DIR="$DATA/outputs/ref_${REF}/cellranger"

OUT="$DATA/outputs/ref_${REF}/seurat"
OUT_map="$OUT/00_mapped"
OUT_qc="$OUT/01_qc"
OUT_filter="$OUT/02_filtered"
OUT_normalize="$OUT/03_normalized"

mkdir -p "$OUT" "$OUT_map" "$OUT_qc" "$OUT_filter" "$OUT_normalize"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Preprocessing for each timepoint in parallel
echo "Preprocessing each sample..."
pids=()
for SMP in "${SMPS[@]}"; do
    mkdir -p "$OUT_qc/$SMP" "$OUT_filter/$SMP"
    srun -N1 -n1 -c 8 --mem=16G "$WORK/scripts/seurat/wrappers/00_preprocessing.sh" "$SMP" &
    pids+=("$!")
done

echo "Waiting for preprocessing jobs to finish..."
for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
        echo "Error: preprocessing failed for at least one sample" >&2
        exit 1
    fi
done

echo "Preprocessing completed✅"


sync
sleep 10


# Integration and clustering all timepoints
echo "Integration and clustering all samples..."
srun -N1 -n1 -c 64 "$WORK/scripts/seurat/wrappers/01_integration_clustering.sh"
echo "Integration and clustering completed✅"


sync
sleep 10


# Sync
echo "Syncing results..."
rsync -av --update --info=progress2 --stats --exclude 'cellranger/' --exclude '*.rds' "$DATA/outputs/" "$WORK/results/"
echo "Syncing completed✅"
