#!/bin/bash
#SBATCH --job-name=map
#SBATCH --partition=courses
#SBATCH --time=01:00:00                 # limit is 24hr
#SBATCH --cpus-per-task=4              	# limit is 128/node
#SBATCH --mem=8G                      	# limit is 515G/node
#SBATCH --chdir=/scratch/dolatabadi.e   # working directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log
#SBATCH --array=0-7

set -euo pipefail

# Arrays
SMPS=("control" "3h" "24h" "72h" "7dpa" "14dpa" "22dpa" "33dpa")

REF="UKY_AmexF1_1_genomic"
SMP=${SMPS[$SLURM_ARRAY_TASK_ID]}

# Paths
DIR="$DATA/outputs/ref_${REF}/cellranger/$SMP/outs/filtered_feature_bc_matrix"
OUT="$DATA/ref_genomes/ref_files"
TSV="$OUT/gene_map.tsv"

cd "$OUT"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Map unknown features
Rscript "$RUN/gene_map.R"
Rscript "$RUN/seurat/0_map_features.R" "$DIR" "$TSV"
