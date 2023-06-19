
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
  setNames(c("12hpi", "24hpi", "4hpi", "72hpi"))

bulk_junc_stats <- map_dfr(bulk_jun_files, read_tsv,
                           col_names = F,
                           show_col_types = F,
                           .id = "timepoint") %>% 
  set_colnames(c("timepoint", "chrom", "start", "end", "junc_name", "read_count", "strand",
                 "thickStart", "thickEnd", "color", "exons", "exon_sizes", "exon_starts"))

bulk_sum_stats <- bulk_junc_stats %>% 
  group_by(timepoint) %>%
  summarize(total_junctions = n(),
            reads10  = sum(read_count > 10),
            reads100  = sum(read_count > 100),
            reads1000  = sum(read_count > 1000),
            total_reads_supporting = sum(read_count),
            mean_read_supporting = total_reads_supporting/total_junctions) %>% 
  mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>% 
  arrange(total_junctions)
  
cov_all <- read_tsv("results/hisat2/coverage/bulk_coverage.txt",
                    comment = "Coverage", show_col_types = FALSE) %>%
  rename("organism" = "#rname") %>%
  filter(organism != "#rname") %>%
  mutate(timepoint = c("4hpi", "12hpi", "24hpi", "72hpi"),
         timepoint = factor(timepoint,
                            levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>%
  map_at(c(2:9), as.numeric) %>%
  as_tibble() %>% 
  select(Timepoint = timepoint, "Total Mapped" = numreads, "Mean Depth" = meandepth)
  
total_seq_reads <- read_table("results/hisat2/coverage/bulk_counts.txt",
                              col_names = T, show_col_types = F)

join_stats <- left_join(cov_all, bulk_sum_stats, by = c("Timepoint" = "timepoint")) %>% 
  bind_cols(total_seq_reads) %>%
  mutate(perc_map = (`Total Mapped`/total_reads) * 100) %>% 
  select(Timepoint, "Total Reads" = total_reads, `Total Mapped`, "Mapped (%)" = perc_map, `Mean Depth`,
         "Splice Junctions" = total_junctions, "Total Junction Reads" = total_reads_supporting,
         reads10, reads100, reads1000) %>%
  map_at(c(2), format, big.mark = ",") %>% as_tibble()

