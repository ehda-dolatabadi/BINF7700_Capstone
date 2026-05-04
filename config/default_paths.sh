#!/bin/bash

# Set the working directory. Modify with absolute path if needed.
export WORK="$(readlink -f .)"

# Set the data directory. Modify if needed.
export DATA="$WORK/data"
export FQ_src="$DATA/Li_dataset"
export FQ_names=("0HPA" "3HPA" "1DPA" "3DPA" "7DPA" "14DPA" "22DPA" "33DPA")
export SMPS=("control" "3h" "24h" "72h" "7dpa" "14dpa" "22dpa" "33dpa")

# Set the directory for reference genome files. Modify if needed.
export REF_FILES="$DATA/ref_files"
export FA_names=("GCF_040938575.1_UKY_AmexF1_1_genomic.fna" "AmexG_v6.0-DD.fa")
export GTF_names=("GCF_040938575.1_UKY_AmexF1_1_genomic.gtf" "AmexT_v47-AmexG_v6.0-DD.noWS.gtf")
export REF_names=("UKY_AmexF1_1_genomic" "AmexT_v47-AmexG_v6_0-DD")

# Set the outputs directory. Modify if needed.
export OUT="$WORK/outputs"

# Set the logs directory. Modify if needed.
export LOG="$WORK/logs"

mkdir -p "$OUT" "$LOG"
