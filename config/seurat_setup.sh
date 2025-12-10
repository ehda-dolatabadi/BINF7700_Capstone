#!/bin/bash
# Automated setup script for Seurat analysis environment

# Step 1: Create conda environment from YAML specification
conda env create -f config/seurat_env.yml

# Step 2: Activate the seurat conda environment
eval "$(conda shell.bash hook)"
conda activate seurat

# Step 3: Restore R packages from renv lockfile
Rscript -e "renv::restore(lockfile = 'config/renv.lock', prompt = FALSE)"
