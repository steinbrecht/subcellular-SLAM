#!/usr/bin/env Rscript

# Load required packages
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# Read input arguments from Snakemake
tsv_files <- snakemake@input$tsv_files
gff3_file <- snakemake@input$gff3_file
output_combined <- snakemake@output$complete
output_gene_level <- snakemake@output$gene_level
output_exon_level <- snakemake@output$exon_level


# Function to process each TSV file
process_tsv_file <- function(file_path) {
  # Read the tsv file
  tsv_data <- read_tsv(file_path)
#   tsv_data <- read_tsv(file_path, col_types = cols(.default = "c"))
  
  # Separate gene and exon IDs
  tsv_data <- tsv_data %>%
    mutate(
    #   gene_id = ifelse(grepl("^ENSMUSG", ID), ID, NA),
      exon_id = ifelse(grepl("^ENSMUSE", ID), ID, "gene")
    )
  
  return(tsv_data)
}

# Read the GFF3 file and extract relationships between gene and exon
gff3_data <- read_tsv(gff3_file, comment = "#", col_names = FALSE) %>%
  filter(grepl("exon", X3)) %>%
  mutate(
    exon_id = sub(".*exon_id=([^;]+);.*", "\\1", X9),
    gene_id_gff = sub(".*gene_id=([^;]+);.*", "\\1", X9)
  ) %>%
  separate(gene_id_gff,c("gene_id_gff","foo"),sep="\\.") %>% 
  separate(exon_id,c("exon_id","bar"),sep="\\.") %>% 
  select(gene_id_gff, exon_id)

print(gff3_data)

# Initialize an empty dataframe to store combined data
combined_data <- NULL

# Process each tsv file
# for (file in readLines(tsv_files)) {
for (file in tsv_files) {
  print(file)

  sample_data <- process_tsv_file(file)
  
  # Merge with gff3 data to get gene-exon relationships
  sample_data <- sample_data %>%
    left_join(gff3_data, by = c("exon_id" = "exon_id")) %>% 
    mutate(gene_id = ifelse(grepl("^ENSMUSG", ID), ID, gene_id_gff),) %>% 
    select(-c(ID, gene_id_gff))
  
  print(sample_data)

  # Combine into a single dataframe
  if (is.null(combined_data)) {
    combined_data <- sample_data
  } else {
    combined_data <- bind_rows(combined_data, sample_data)
  }
}


# put the columns in good order
combined_data <- combined_data %>% 
    select(c(sample, gene_id, exon_id), everything())


# ---- COMBINED -------------------
# Write the combined data to the output file
write_tsv(combined_data, output_combined)


# ---- ONLY GENE LEVEL -------------------
gene_data <- combined_data %>% 
    filter(grepl("gene",exon_id)) %>% 
    select(-exon_id)

write_tsv(gene_data, output_gene_level)


# ---- ONLY EXON LEVEL -------------------
# this takes the intronic T2C data for the gene-level entries 
# and the exonic T2C data for the exonic entries
exon_data <- combined_data %>% 
    mutate(
        T2C_over_T_normalized = ifelse(grepl("gene",exon_id), T2C_over_T_intron, T2C_over_T_exon),
        T_count = ifelse(grepl("gene",exon_id), T_intron, T_exon),
    ) %>% 
    select(-c(T2C_over_T_intron, T2C_over_T_exon, T_intron, T_exon, T2C_intron, T2C_exon))

write_tsv(exon_data, output_exon_level)