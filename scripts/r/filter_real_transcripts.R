#!/usr/bin/env Rscript


### FILTER REAL TRANSCRIPTS FROM PREDICTED ORFS AND REMOVE DUPLICATES############

library(glue)
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

# function to filter any given time-point
sub_timepoint <- function(timepoint){
  out <- merged_gtf %>%
    filter(str_detect(X9, "ref_gene_name", negate = T)) %>% 
    filter(str_detect(X9, paste0('gene_id\\s"', timepoint)))
  
  # save filtered results
  write.table(out, file = paste0("results/stringtie/transcripts_merged_", timepoint, "hrs.gtf"),
              col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}

sub_timepoint(4)
sub_timepoint(12)
sub_timepoint(24)
sub_timepoint(72)
