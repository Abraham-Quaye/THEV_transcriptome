#!/usr/bin/env Rscript

library(SGSeq)
library(tidyverse)

## make list of samples with corresponding bam files to extract information
# sample <- tibble(sample_name = c(paste0("12hrs", c("S1", "S3")),
#                                  paste0(c("24hrs"), c("S1", "S2", "S3")),
#                                  paste0(c("4hrs"), c("S1", "S2", "S3")),
#                                  paste0(c("72hrs"), c("S1", "S2", "S3"))),
#                  file_bam = list.files("results/hisat2", "^thev_sorted_.+bam$",
#                                        full.names = T)
#                  )

# ## extract information from bam files. This is a slow function so I'll save the data
# # as a tsv and load it henceforth
# sample_info <- getBamInfo(sample, cores = 8)

# write.table(sample_info, "results/sgseq/bam_extracted_sampleData.tsv",
#             col.names = T, row.names = F, quote = F, sep = "\t")

### ====================================
# load in sample/bam information

bam_info <- read_tsv("results/sgseq/bam_extracted_sampleData.tsv",
  show_col_types = F) 

bam_info %>%
  select(read_length) %>%
  slice(1) %>%
  pull()

