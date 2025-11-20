#!/bin/bash

set -euo pipefail

# ==================== VARIABLES ====================
DATA="/courses/BINF7700.202610/students/dolatabadi.e"
SCRATCH="/scratch/dolatabadi.e"
WORK="$HOME/BINF7700_Capstone"
RUN="$WORK/scripts"
SLURM="$WORK/scripts/seurat/slurm"
TSV="$DATA/ref_genomes/ref_files/loc_map.tsv"
OUT="$DATA/outputs/seurat"

main_ID="main"

cluster_ID="cluster14"
cluster="14"

cells_ID="epith-eryth"

# ==================== EXPORTS ====================
export DATA SCRATCH WORK RESULT RUN SLURM TSV OUT
export main_ID cluster_ID cells_ID IDENT

mkdir -p "$OUT" "$SCRATCH/logs"

# ==================== PIPELINE ====================
# Main Pipeline
# Step 1: Preprocessing
echo "Submitting preprocessing..."
#JOB1=$(sbatch --parsable --export=ALL,ID=$main_ID $SLURM/01_preprocess.sbatch)
#echo "  Job ID: $JOB1"

# Step 2: Integration
echo "Submitting integration..."
#JOB2=$(sbatch --parsable --export=ALL,ID=$main_ID --dependency=afterok:$JOB1 $SLURM/02_integrate.sbatch)
#echo "  Job ID: $JOB2"

# Step 3: Clustering
echo "Submitting clustering..."
#JOB3=$(sbatch --parsable --export=ALL,ID=$main_ID --dependency=afterok:$JOB2 $SLURM/03_cluster.sbatch)
#echo "  Job ID: $JOB3"



# Subclustering
echo "Submitting subclustering..."
JOB4=$(sbatch --parsable --export=ALL,ID=$cluster_ID,IDENT=$cluster $SLURM/04_subcluster.sbatch) #--dependency=afterok:$JOB3 $SLURM/04_subcluster.sbatch)
echo "  Job ID: $JOB4"

echo "Submitting clustering for subset..."
JOB5=$(sbatch --parsable --export=ALL,ID=$cluster_ID --dependency=afterok:$JOB4 $SLURM/03_cluster.sbatch)
echo "  Job ID: $JOB5"



# Removing abundant cells
echo "Submitting removal..."
#JOB6=$(sbatch --parsable --export=ALL,ID=$cells_ID $SLURM/05_subcells.sbatch) #--dependency=afterok:$JOB2 $SLURM/05_subcells.sbatch)
#echo "  Job ID: $JOB6"

echo "Submitting clustering for subset..."
#JOB7=$(sbatch --parsable --export=ALL,ID=$cells_ID --dependency=afterok:$JOB6 $SLURM/03_cluster.sbatch)
#echo "  Job ID: $JOB7"
echo "Pipeline submitted successfully!"
