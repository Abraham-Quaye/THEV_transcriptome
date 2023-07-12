#!/usr/bin/env Rscript


### FILTER REAL TRANSCRIPTS FROM PREDICTED ORFS AND REMOVE DUPLICATES############

library(tidyverse)

# Merged all StringTie gff files for all timepoints -> all_merged.gtf --------

# Remove all predicted genes from all_merged.gff file to view only 
# real transcripts and exons in IGV and help plan the primers 
# for wet-lab transcript validation
##### Also remove duplicates
# load gtf file ------------
merged_gtf <- read_tsv("results/stringtie/all_merged.gtf",
                       comment = "#",
                       col_names = FALSE,
                       show_col_types = FALSE) 

duplicates <- 'gene_id "24S2.11"|gene_id "24S1.11"|gene_id "72S3.10"|transcript_id "24S1.2.1"|transcript_id "24S3.2.5"|gene_id "12S1.2"|transcript_id "24S1.1.1"|gene_id "24S2.1"|transcript_id "24S2.2.5"|transcript_id "72S1.8.1"|transcript_id "72S2.8.1"|transcript_id "4S2.3.1"|transcript_id "72S1.10.1"|transcript_id "24S1.12.1"|transcript_id "72S2.9.1"|transcript_id "72S1.1.1"|transcript_id "72S1.2.2"|transcript_id "72S3.1.1"|transcript_id "24S3.1.1"|transcript_id "24S1.1.2"|transcript_id "12S1.1.1"|transcript_id "24S2.10.1"|transcript_id "12S1.10.2"|transcript_id "24S2.2.1"|transcript_id "4S1.4.1"'

filter_real_transcripts <- merged_gtf %>% 
  filter(!str_detect(X9,"ref_gene_name")) %>% 
  filter(!str_detect(X9, duplicates)) %>%
  arrange(X4)

# save modified .gtf
write.table(filter_real_transcripts, "results/stringtie/all_real_transcripts_merged.gtf",
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


