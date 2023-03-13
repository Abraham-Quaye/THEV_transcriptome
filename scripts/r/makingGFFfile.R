#!/usr/bin/env Rscript

library(tidyverse)

# bedfile <- read_tsv("raw_files/annotations/THEVannotated_genesOnly.txt",  col_names = F, col_types = cols(X11 = "c"))
# write.table(bedfile, "THEVannotated_genesOnly.bed", col.names = F, row.names = F, sep = "\t", quote = F)

## agat used to convert .bed to .gff

# .gff contains 3' and 5' UTRs that have the startsite downstream of the stopsite -> errors.
# I'm removing all rows with such information

# load gff file
gfffile <- read_tsv("raw_files/annotations/thev_predicted_genes.gff", col_names = F, skip = 1) %>%
    filter(!X3 %in% c("three_prime_UTR","five_prime_UTR"))

# replace old .gff with modified .gff
write.table(gfffile, "raw_files/annotations/thev_predicted_genes.gff", col.names = F, row.names = F, sep = "\t", quote = F)



