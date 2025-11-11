#!/usr/bin/env Rscript

# Args
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
tsv <- args[2]

input <- file.path(dir, "features.tsv.gz")

df <- read.table(tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")

features <- read.table(gzfile(input), 
                       sep = "\t", 
                       header = FALSE, 
                       stringsAsFactors = FALSE)

colnames(features) <- c("id", "name", "type")

gene_map <- setNames(df$symbol, df$gene_id)
features$name <- ifelse(features$name %in% names(gene_map),
                        gene_map[features$name],
                        features$name)

write.table(features,
            gzfile(input),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
