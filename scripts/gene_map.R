#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

dir <- "/courses/BINF7700.202610/students/dolatabadi.e/ref_genomes/ref_files"
gtf <- file.path(dir, "GCF_040938575.1_UKY_AmexF1_1_genomic.gtf")
hgnc <- file.path(dir, "hgnc_complete_set.tsv")
manual <- file.path(dir, "manual_map.tsv")
out <- file.path(dir, "gene_map.tsv")

# read HGNC
HGNC <- read_tsv(
  hgnc,
  quote = "",
  comment = "",
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
)

# read and filter genes from gtf
genes <- read.table(gtf, sep = "\t") %>% 
  filter(V3 == "gene") %>%
  filter(str_detect(V9, "description")) %>% 
  select(V9) %>% 
  mutate(gene_id = sapply(str_split(V9, ";"), function(x) x[1])) %>% 
  mutate(gene_id = str_remove(gene_id, "gene_id ")) %>% 
  mutate(description = sapply(str_split(V9, ";"), function(x) x[4])) %>% 
  mutate(description = str_remove(description, " description ")) %>% 
  select(gene_id, description) %>% 
  arrange(gene_id)

# read manuallly mapped genes
manual <- read.table(manual, header = TRUE, sep = "\t")

# filter characterized loc genes 
loc_genes <- genes %>%
  filter(str_starts(gene_id, "LOC")) %>% 
  filter(!str_detect(description, "uncharacterized "))

# define filler words
filler_words <- c("associated", 
                  "containing", 
                  "family", 
                  "homolog", 
                  "isoform", 
                  "ligand",
                  "like",
                  "member",
                  "membrane",
                  "precursor", 
                  "probable", 
                  "pseudogene",
                  "putative", 
                  "related",  
                  "subunit",
                  "type")

# normalize function
normalize <- function(text) {
  text %>%
    tolower() %>%
    str_replace_all("\\bmitochondrial-like\\b", "") %>% 
    str_replace_all("[^a-z0-9\\s]", " ") %>% 
    str_replace_all("\\bmitochondrial pseudogene\\b", "") %>% 
    str_replace_all(paste0("\\b(", paste(filler_words, collapse = "|"), ")\\b"), "") %>% 
    str_replace_all("\\s+", " ") %>% 
    str_trim()
}

# normalize loc genes
loc_norm <- loc_genes %>%
  arrange(description) %>% 
  distinct(description) %>% 
  mutate(description_norm = normalize(description))

# normalize defined genes in gtf
genes_norm <- genes %>% 
  filter(!str_starts(gene_id, "LOC")) %>%
  mutate(gene_id = 
           if_else(
             duplicated(select(., description)) | duplicated(select(., description), fromLast = TRUE),
             str_sub(gene_id, 1, pmin(9, str_length(gene_id))),
             gene_id)
         ) %>% 
  distinct(gene_id, description, .keep_all = TRUE) %>% 
  mutate(description_norm = normalize(description))

# match loc genes with defined genes
matched1 <- inner_join(loc_norm, genes_norm, by = c("description_norm" = "description_norm"), keep = TRUE) %>% 
  select(description.x, gene_id) %>%
  distinct(description.x, .keep_all = TRUE) %>% 
  rename(symbol = gene_id)

# exclude define genes
loc_norm <- loc_norm %>%
  filter(!description %in% matched1$description.x)

# match to loc ids
matched1 <- inner_join(loc_genes, matched1, by = c("description" = "description.x"))

# normalize HGNC
HGNC_norm <- HGNC %>%
  select(symbol, name, alias_name, prev_name) %>% 
  arrange(symbol) %>%
  mutate(
    alias_name = str_split(coalesce(alias_name, ""), "\\|"),
    prev_name  = str_split(coalesce(prev_name, ""), "\\|")
  ) %>% 
  rowwise() %>%
  reframe(symbol = symbol, name = c(name, alias_name, prev_name)) %>%
  filter(name != "") %>%
  mutate(name_norm = normalize(name)) %>% 
  distinct(name, .keep_all = TRUE)

# match loc genes with HGNC
matched2 <- inner_join(loc_norm, HGNC_norm, by = c("description_norm" = "name_norm"), keep = TRUE) %>% 
  select(description, symbol) %>%
  distinct(description, .keep_all = TRUE)

# exclude defined genes
loc_norm <- loc_norm %>%
  filter(!description %in% matched2$description)

# match to loc ids
matched2 <- inner_join(loc_genes, matched2, by = c("description" = "description"))

# normalize manually mapped genes
manual_norm <- manual %>% 
  mutate(description_norm = normalize(description))

matched3 <- inner_join(loc_norm, manual_norm, by = c("description_norm" = "description_norm"), keep = TRUE) %>% 
  select(description.x, symbol) %>%
  distinct(description.x, .keep_all = TRUE)

# exclude defined genes
loc_norm <- loc_norm %>%
  filter(!description %in% matched3$description.x)

# match to loc ids
matched3 <- inner_join(loc_genes, matched3, by = c("description" = "description.x"))

# combine all
matched <- manual %>% 
  bind_rows(matched1) %>% 
  bind_rows(matched2) %>% 
  bind_rows(matched3) %>% 
  distinct(gene_id, .keep_all = TRUE) %>% 
  arrange(gene_id)

# write to tsv
write.table(matched, file=out, sep = "\t", quote = FALSE, row.names = FALSE)
