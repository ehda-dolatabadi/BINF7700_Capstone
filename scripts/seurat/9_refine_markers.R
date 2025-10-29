#!/usr/bin/env Rscript
# Purpose: Refine marker list by filtering and keeping top informative genes per cluster
# Usage: Rscript 9_refine_markers.R <id> <outdir> <input_markers>

suppressPackageStartupMessages({
  library(dplyr)
})

set.seed(777)

args <- commandArgs(trailingOnly = TRUE)
id     <- args[1]
outdir <- args[2]
input  <- args[3]

markers <- read.table(input, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Filter
markers <- markers %>%
  filter(
	p_val_adj < 0.05,
	avg_log2FC >= 0.5,
	(pct.1 - pct.2) >= 0.2
  )

# Rank and pick top 10 per cluster
top <- markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), desc(pct.1 - pct.2)) %>%
  slice_head(n = 10) %>%
  ungroup()

outfile <- file.path(outdir, paste0(id, "_markers_refined.tsv"))
write.table(top, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
