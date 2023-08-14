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

duplicates <- 'gene_id "24S2.11"|gene_id "24S1.11"|gene_id "72S3.10"|transcript_id "24S3.2.5"|gene_id "12S1.2"|transcript_id "24S1.1.1"|gene_id "24S2.1"|transcript_id "24S2.2.5"|transcript_id "72S1.8.1"|transcript_id "72S2.8.1"|transcript_id "4S2.3.1"|transcript_id "72S1.10.1"|transcript_id "24S1.12.1"|transcript_id "72S2.9.1"|transcript_id "72S1.1.1"|transcript_id "72S3.1.1"|transcript_id "24S3.1.1"|transcript_id "24S1.1.2"|transcript_id "12S1.1.1"|transcript_id "24S2.10.1"|transcript_id "12S1.10.2"|transcript_id "24S2.2.1"|transcript_id "4S1.4.1"|transcript_id "24S2.2.2"|transcript_id "24S3.2.1"|transcript_id "72S1.12.1"|transcript_id "24S1.2.4"|transcript_id "12S1.11.3"|transcript_id "24S3.11.1"|transcript_id "72S1.13.1"|transcript_id "72S1.13.2"|transcript_id "24S1.14.1"|transcript_id "12S1.4.2"'

filter_real_transcripts <- merged_gtf %>% 
  filter(str_detect(X9, "reference_id", negate = T)) %>%
  filter(!str_detect(X9, duplicates)) %>%
  group_by(X7) %>% 
  arrange(X4, .by_group = T)

# save modified .gtf of all unique transcripts across all time-points
write.table(filter_real_transcripts, "results/stringtie/all_real_transcripts_merged.gtf",
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


# function to filter any given time-point
sub_timepoint <- function(timepoint){
  out <- merged_gtf %>%
    filter(str_detect(X9, "reference_id", negate = T)) %>%
    filter(str_detect(X9, paste0('gene_id\\s"', timepoint)))
  
  # save filtered results
  write.table(out, file = paste0("results/stringtie/transcripts_merged_", timepoint, "hrs.gtf"),
              col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}

sub_timepoint(4)
sub_timepoint(12)
sub_timepoint(24)
sub_timepoint(72)
