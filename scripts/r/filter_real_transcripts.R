#!/usr/bin/env Rscript


### FILTER REAL TRANSCRIPTS FROM PREDICTED ORFS ############

library(tidyverse)

# Merged all StringTie gff files for all timepoints -> all_merged.gtf --------

# Remove all predicted genes from all_merged.gff file to view only 
# real transcripts and exons in IGV and help plan the primers 
# for wet-lab transcript validation

# load gtf file ------------
merged_gtf <- read_tsv("results/stringtie/all_merged.gtf",
                       comment = "#",
                       col_names = FALSE,
                       show_col_types = FALSE) 

filter_real_transcripts <- merged_gff %>% 
  filter(!str_detect(X9,"ref_gene_name"))

# replace old merged .gff with modified .gff
write.table(filter_real_transcripts, "results/stringtie/all_real_transcripts_merged.gtf",
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
