#!/bin/bash
#SBATCH --job-name=download2
#SBATCH --partition=courses
#SBATCH --time=24:00:00			# limit 24hr
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --chdir=/scratch/dolatabadi.e   # working directory
#SBATCH --output=logs/%x_%j.out         # stdout log
#SBATCH --error=logs/%x_%j.err          # stderr log

set -euo pipefail

#rsync -av --update --info=progress2 --stats "$SCRATCH/outputs/" "$DATA/outputs/"

#rsync -av --update --info=progress2 --stats "$SCRATCH/ref_genomes/cellranger_mkref/" "$DATA/ref_genomes/cellranger_mkref/"

#rsync -av --ignore-existing --info=progress2 --stats "$DATA/3HPA" "$SCRATCH/Li_dataset/"
#rsync -av --partial --append-verify --info=progress2 --stats "$DATA/3HPA" "$SCRATCH/Li_dataset/"

#mkdir -p $DATA/downloads/
#cd $DATA/downloads/
#wget -c ftp://ftp2.cngb.org/pub/CNSA/data2/CNP0000706/CNS0130458/CNX0115114/CNR0137411/3h_S0_L001_R2_001.fastq.gz
#curl -C - -O --retry 100 --retry-delay 30 --speed-limit 1000 --speed-time 60 \
  #ftp://ftp2.cngb.org/pub/CNSA/data2/CNP0000706/CNS0130458/CNX0115114/CNR0137411/3h_S0_L001_R2_001.fastq.gz

rsync -av --update --info=progress2 --stats "$SCRATCH/outputs/cellranger/ref_AmexT_v47-AmexG_v6_0-DD/14dpa" "$DATA/outputs/cellranger/ref_AmexT_v47-AmexG_v6_0-DD/"
