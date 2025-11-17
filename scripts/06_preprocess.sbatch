#!/bin/bash
#SBATCH --job-name=preprocess
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
REFS=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")

SMPS=(		"control"	"3h"	"24h"	"72h"	"7dpa"	"14dpa"	"22dpa"	"33dpa")
min_counts=(	"1000" 		"1000" 	"1000" 	"1000" 	"1000" 	"1000" 	"1000" 	"1000")
#max_counts=(	"20000"		"12500"	"10000"	"20000"	"30000"	"18000"	"30000"	"25000")
max_counts=(	"30000"         "30000" "30000" "30000" "30000" "30000" "30000" "30000")
min_features=(	"500"		"500"	"500"	"500"	"500"	"500"	"500"	"500")
max_features=(	"5000"		"5000"	"5000"	"5000"	"5000"	"5000"	"5000"	"5000")
max_mt=(	"10"		"10"	"10"	"10"	"10"	"10"	"10"	"10")
max_ribo=(	"30"		"30"	"30"	"30"	"30"	"30"	"30"	"30")

# Derive indices
ref_idx=$((SLURM_ARRAY_TASK_ID / ${#SMPS[@]}))
smp_idx=$((SLURM_ARRAY_TASK_ID % ${#SMPS[@]}))

REF=${REFS[$ref_idx]}
SMP=${SMPS[$smp_idx]}

p1=${min_counts[$smp_idx]}
p2=${max_counts[$smp_idx]}
p3=${min_features[$smp_idx]}
p4=${max_features[$smp_idx]}
p5=${max_mt[$smp_idx]}
p6=${max_ribo[$smp_idx]}

# Paths
DIR="$DATA/outputs/ref_${REF}/seurat/00_mapped/${SMP}_mapped.rds"
OUT_filter="$DATA/outputs/ref_${REF}/seurat/02_filtered/$SMP"
OUT_normalize="$DATA/outputs/ref_${REF}/seurat/03_normalized"

mkdir -p "$OUT_filter"
mkdir -p "$OUT_normalize"

# Environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Run filter
Rscript "$RUN/seurat/02_filter.R" "$SMP" "$OUT_filter" "$DIR"\
  "$p1" "$p2" "$p3" "$p4" "$p5" "$p6"

# Run normalize
Rscript "$RUN/seurat/03_normalize.R" "$SMP" "$OUT_normalize" "$OUT_filter/${SMP}_filtered.rds" "${SLURM_CPUS_PER_TASK}"

# Sync
rsync -av --update --info=progress2 --stats --exclude 'cellranger/' --exclude '*.rds' "$DATA/outputs/" "$WORK/results/"
