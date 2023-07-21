
library(magrittr)
library(glue)
library(tidyverse)


single_jun_files <- list.files("results/hisat2", full.names = T,
                               pattern = "\\w+\\.bed") %>% 
  setNames(c("t12s1","t12s3", "t24s1", "t24s2", "t24s3", "t4s1", "t4s2", "t4s3", "t72s1", "t72s2", "t72s3"))

single_junc_stats <- map_dfr(single_jun_files, read_tsv,
                             col_names = F,
                             show_col_types = F,
                             .id = "timepoint") %>% 
  set_colnames(c("timepoint", "chrom", "start", "end", "junc_name", "read_count", "strand",
                 "thickStart", "thickEnd", "color", "exons", "exon_sizes", "exon_starts")) %>% 
  arrange(start, end)

single_junclist <- single_junc_stats %>% group_by(start, end) %>% summarise(repeated = n(),
                                                                            .groups = "drop")
 
single_sum_stats <- single_junc_stats %>% 
  mutate(timepoint = case_when(timepoint %in% c("t72s1", "t72s2", "t72s3") ~ "72hpi",
                               timepoint %in% c("t24s1", "t24s2", "t24s3") ~ "24hpi",
                               timepoint %in% c("t12s1", "t12s3") ~ "12hpi",
                               timepoint %in% c("t4s1", "t4s2", "t4s3") ~ "4hpi",
                               .default = timepoint)) %>% 
  group_by(timepoint) %>%
  summarize(total_junctions = n(),
            reads10  = sum(read_count > 10),
            reads100  = sum(read_count > 100),
            reads1000  = sum(read_count > 1000),
            total_reads_supporting = sum(read_count),
            mean_read_supporting = total_reads_supporting/total_junctions,
            organism = "thev") %>%
  # mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>% 
  arrange(total_junctions)

# ============================BULK=============================================
bulk_jun_files <- list.files("results/hisat2/bulk", full.names = T,
                             pattern = "\\w+\\.bed") %>% 
  setNames(c("12hpi", "24hpi", "4hpi", "72hpi"))

bulk_junc_stats <- map_dfr(bulk_jun_files, read_tsv,
                           col_names = F,
                           show_col_types = F,
                           .id = "timepoint") %>% 
  set_colnames(c("timepoint", "chrom", "start", "end", "junc_name", "read_count", "strand",
                 "thickStart", "thickEnd", "color", "exons", "exon_sizes", "exon_starts")) %>% 
  arrange(start, end)

# bulk_junclist <- bulk_junc_stats %>% group_by(start, end) %>% summarise(repeated = n(),
#                                                                         .groups = "drop")

bulk_sum_stats <- bulk_junc_stats %>% 
  group_by(timepoint) %>%
  summarize(total_junctions = n(),
            reads10  = sum(read_count > 10),
            reads100  = sum(read_count > 100),
            reads1000  = sum(read_count > 1000),
            total_reads_supporting = sum(read_count),
            mean_read_supporting = total_reads_supporting/total_junctions,
            organism = "thev") %>% 
  mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>% 
  arrange(total_junctions)

thev_totals <- tibble(sum_junc = sum(bulk_sum_stats$total_junctions), sum_reads10 = sum(bulk_sum_stats$reads10),
            sum_reads100 = sum(bulk_sum_stats$reads100), sum_reads1000 = sum(bulk_sum_stats$reads1000),
            sum_junc_reads = sum(bulk_sum_stats$total_reads_supporting)) %>% 
  map_at(c(1:5), format, big.mark = ",", drop0trailing = T) %>% as_tibble()

bulk_sum_stats <- bulk_sum_stats %>%
  map_at(c(2:6), format, big.mark = ",", drop0trailing = T, nsmall = 2) %>%
  as_tibble()

total_reads <- read_table("results/hisat2/coverage/bulk_counts.txt",
                          col_names = T, show_col_types = F)

### find all the coverage files
all_cov_files <- list.files("results/hisat2/coverage",
                            pattern = "host_thev_cov\\d+.txt",
                            full.names = TRUE) %>% 
  setNames(c("12hpi", "24hpi", "4hpi", "72hpi"))

# read depth_files as one master tibble
all_covs <- map_dfr(all_cov_files, read_tsv,
                    show_col_types = FALSE,
                    comment = "Coverage",
                    .id = "timepoint") %>%
  rename("organism" = "#rname") %>% 
  mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>% 
  map_at(c(3:10), as.numeric) %>% 
  as_tibble() %>%
  mutate(organism = ifelse(organism == "AY849321.1", "thev", "m.gallopavo")) %>% 
  group_by(organism, timepoint) %>%
  reframe(total_mapped = sum(numreads),
          mean_depth = mean(meandepth), mean_cov = mean(coverage)) %>% 
  cbind(total_reads) %>% 
  mutate(perc_mapped = round((total_mapped/total_reads)*100,4)) %>%
  map_at(c(3:7), format, big.mark = ",", digits = 3, drop0trailing = T) %>%
  as_tibble()

join_stats <- left_join(all_covs, bulk_sum_stats, by = c("timepoint", "organism")) %>%
  select(organism, timepoint, total_reads, total_mapped, perc_mapped, mean_depth, mean_cov,
         total_junctions, total_reads_supporting, reads10, reads100, reads1000)


tab1 <- tibble(Metric = c("Total reads", "Mapped\n(Host)", "Mapped\n(THEV)", "Splice junctions", "Junction coverage\n>= 1 read","Junction coverage\n>= 10 reads","Junction coverage\n>= 100 reads", "Junction coverage\n>= 1000 reads"),
               "4h.p.i" = c(all_covs$total_reads[1], glue("{all_covs$total_mapped[1]} ({all_covs$perc_mapped[1]}%)"), glue("{all_covs$total_mapped[5]} ({all_covs$perc_mapped[5]}%)"), join_stats$total_junctions[5], join_stats$total_reads_supporting[5], join_stats$reads10[5], join_stats$reads100[5], join_stats$reads1000[5]),
               "12h.p.i" = c(all_covs$total_reads[2], glue("{all_covs$total_mapped[2]} ({all_covs$perc_mapped[2]}%)"), glue("{all_covs$total_mapped[6]} ({all_covs$perc_mapped[6]}%)"), join_stats$total_junctions[6], join_stats$total_reads_supporting[6], join_stats$reads10[6], join_stats$reads100[6], join_stats$reads1000[6]),
               "24h.p.i" = c(all_covs$total_reads[3], glue("{all_covs$total_mapped[3]} ({all_covs$perc_mapped[3]}%)"), glue("{all_covs$total_mapped[7]} ({all_covs$perc_mapped[7]}%)"), join_stats$total_junctions[7], join_stats$total_reads_supporting[7], join_stats$reads10[7], join_stats$reads100[7], join_stats$reads1000[7]),
               "72h.p.i" = c(all_covs$total_reads[4], glue("{all_covs$total_mapped[4]} ({all_covs$perc_mapped[4]}%)"), glue("{all_covs$total_mapped[8]} ({all_covs$perc_mapped[8]}%)"), join_stats$total_junctions[8], join_stats$total_reads_supporting[8], join_stats$reads10[8], join_stats$reads100[8], join_stats$reads1000[8]),
               Total = c(sum(as.numeric(all_covs$total_reads[1:4])), sum(as.numeric(all_covs$total_mapped[1:4])), sum(as.numeric(all_covs$total_mapped[5:8])), thev_totals$sum_junc, thev_totals$sum_junc_reads, thev_totals$sum_reads10, thev_totals$sum_reads100, thev_totals$sum_reads1000)) 

