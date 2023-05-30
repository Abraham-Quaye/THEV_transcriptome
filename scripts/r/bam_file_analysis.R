
library(magrittr)
library(tidyverse)


single_jun_files <- list.files("results/hisat2", full.names = T,
                               pattern = "\\w+\\.bed") %>% 
  setNames(c("t12s1","t12s3", "t24s1", "t24s2", "t24s3", "t4s1", "t4s2", "t4s3", "t72s1", "t72s2", "t72s3"))

single_junc_stats <- map_dfr(single_jun_files, read_tsv,
                             col_names = F,
                             show_col_types = F,
                             .id = "timepoint") %>% 
  set_colnames(c("timepoint", "chrom", "start", "end", "junc_name", "read_count", "strand",
                 "thickStart", "thickEnd", "color", "exons", "exon_sizes", "exon_starts"))

bulk_jun_files <- list.files("results/hisat2/bulk", full.names = T,
                             pattern = "\\w+\\.bed") %>% 
  setNames(c("t12hrs", "t24hrs", "t4hrs", "t72hrs"))

bulk_junc_stats <- map_dfr(bulk_jun_files, read_tsv,
                           col_names = F,
                           show_col_types = F,
                           .id = "timepoint") %>% 
  set_colnames(c("timepoint", "chrom", "start", "end", "junc_name", "read_count", "strand",
                 "thickStart", "thickEnd", "color", "exons", "exon_sizes", "exon_starts"))

bulk_sum_stats <- bulk_junc_stats %>% 
  group_by(timepoint) %>%
  summarize(total_junctions = n(),
            total_reads_supporting = sum(read_count),
            mean_read_supporting = total_reads_supporting/total_junctions)
  
stats_per_timepoint
  
  
  
