#!/usr/bin/env Rscript
# Script: install_packages.R
# Purpose: Install CellChat and GitHub-only dependencies into the conda environment
# Description: Run once before the per-timepoint array job. Idempotent — skips
#              packages already installed. Per-timepoint scripts only do library()
#              calls; they assume the environment is already set up.
# Usage: Rscript install_packages.R

suppressPackageStartupMessages({
  library(remotes)
})

if (!requireNamespace("NMF", quietly = TRUE))
    remotes::install_version("NMF", version = "0.28.0",
                             repos = "https://cloud.r-project.org",
                             upgrade = "never")

if (!requireNamespace("CellChat", quietly = TRUE))
    remotes::install_github("jinworks/CellChat@v2.1.2", upgrade = "never")

if (!requireNamespace("presto", quietly = TRUE)) {
    Sys.setenv("PKG_CXXFLAGS" = "-std=gnu++14")
    remotes::install_github("immunogenomics/presto@1.0.0", upgrade = "never")
}

