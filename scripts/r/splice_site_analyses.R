#!/usr/bin/env Rscript

library(magrittr)
library(scales)
library(tidyverse)
library(rtracklayer)
library(glue)
library(ggsci)

# load in data from from full transcriptome to access the splice sites
all_spliced_gtf <- import("results/gffcompare/gffcomp_alltimes.combined.gtf") %>%
  as_tibble() %>%
  select(chr = seqnames, start, end, orf_size = width, strand, type, gene_id,
         transcript_id, exon_number) %>% 
  group_by(transcript_id, strand) %>%
  reframe(start = list(start),
          end = list(end)) %>%
  unnest_wider(c(start, end), names_sep = "_") %>%
  mutate(region = case_when(transcript_id %in% paste0("TCONS_000000", c('15', '16', '04', '05')) ~ "E1",
                            transcript_id %in% paste0("TCONS_000000", c('03', '17', '18', '28')) ~ "E2",
                            transcript_id %in% paste0("TCONS_000000", c('26', '27', '12', '10', '11', '02')) ~ "E3",
                            transcript_id == "TCONS_00000014" ~ "E4",
                            transcript_id == "TCONS_00000013" ~ "IM",
                            TRUE ~ "MLP")) %>% 
  arrange(strand, region) %>%
  mutate(j1_ss = end_2, j1_ts = start_3,
         j2_ss = end_3, j2_ts = start_4,
         j3_ss = end_4, j3_ts = start_5,
         j4_ss = end_5, j4_ts = start_6,
         j5_ss = end_6, j5_ts = start_7,
         j6_ss = end_7, j6_ts = start_8
         ) %>% 
  select(-ends_with("_1"), -starts_with("start"), -starts_with("end"))

# extract junction pairs (start and end) to make long format
get_junc_pair <- function(j_ss, j_ts){
  return(data.frame(transcript_id = all_spliced_gtf$transcript_id,
             strand = all_spliced_gtf$strand,
             region = all_spliced_gtf$region,
             junc_ss = all_spliced_gtf %>% select(all_of(j_ss)),
             junc_ts = all_spliced_gtf %>% select(all_of(j_ts))
             ) %>%
           set_colnames(c("transcript_id", "strand", "region",
                                  "junc_ss", "junc_ts")) %>%
           as_tibble()
  )
}


tot_unique_spliced_juncs <- rbind(get_junc_pair("j1_ss", "j1_ts"), get_junc_pair("j2_ss", "j2_ts"),
                                  get_junc_pair("j3_ss", "j3_ts"), get_junc_pair("j4_ss", "j4_ts"),
                                  get_junc_pair("j5_ss", "j5_ts"), get_junc_pair("j6_ss", "j6_ts")) %>%
  drop_na() %>% 
  distinct(junc_ss, junc_ts, .keep_all = T) %>%
  mutate(junc_name = glue("junc_{c(1:nrow(.))}"))


# load in counted splice sites from regtools output to quantify 
# the splices sites of any specified region 
# -----------------
single_jun_files <- list.files("results/hisat2", full.names = T,
                               pattern = "\\w+\\.bed") %>%
  setNames(c("t12s1","t12s3", "t24s1", "t24s2", "t24s3","t4s1",
             "t4s2", "t4s3", "t72s1", "t72s2", "t72s3"))

single_junc_stats <- map_dfr(single_jun_files, read_tsv,
                             col_names = F,
                             show_col_types = F,
                             col_types = c("ciiciciicicc"),
                             .id = "timepoint") %>%
  set_colnames(c("timepoint", "chrom", "start", "end", "junc_name",
                 "read_count", "strand", "thickStart", "thickEnd",
                 "color", "exons", "exon_sizes", "exon_starts")) %>%
  arrange(start, end) %>%
  distinct(start, end, .keep_all = T) %>%
  mutate(timepoint = case_when(timepoint %in% glue("t4s{c(1,2,3)}") ~ "4hpi",
                               timepoint %in% glue("t12s{c(1,3)}") ~ "12hpi",
                               timepoint %in% glue("t24s{c(1,2,3)}") ~ "24hpi",
                               timepoint %in% glue("t72s{c(1,2,3)}") ~ "72hpi",
                               TRUE ~ timepoint)) %>% 
  separate(exon_sizes, c("exon1size", "exon2size"), sep = ",", convert = T) %>%
  mutate(junc_ss = start + exon1size,
         junc_ts = (end - exon2size) + 1)


# -----------------
bulk_jun_files <- list.files("results/hisat2/bulk", full.names = T,
                             pattern = "\\w+\\.bed") %>%
  setNames(c("12hpi", "24hpi", "4hpi", "72hpi"))

bulk_junc_stats <- map_dfr(bulk_jun_files, read_tsv,
                           col_names = F,
                           col_types = c("ciiciciicicc"),
                           show_col_types = F,
                           .id = "timepoint") %>%
  set_colnames(c("timepoint", "chrom", "start", "end", "junc_name",
                 "read_count", "strand", "thickStart", "thickEnd",
                 "color", "exons", "exon_sizes", "exon_starts")) %>%
  distinct(start, end, .keep_all = T) %>%
  arrange(start, end) %>%
  distinct(start, end, .keep_all = T) %>%
  mutate(timepoint = case_when(timepoint %in% glue("t4s{c(1,2,3)}") ~ "4hpi",
                               timepoint %in% glue("t12s{c(1,3)}") ~ "12hpi",
                               timepoint %in% glue("t24s{c(1,2,3)}") ~ "24hpi",
                               timepoint %in% glue("t72s{c(1,2,3)}") ~ "72hpi",
                               TRUE ~ timepoint)) %>% 
  separate(exon_sizes, c("exon1size", "exon2size"), sep = ",", convert = T) %>%
  mutate(junc_ss = start + exon1size,
         junc_ts = (end - exon2size) + 1)
  


# filter out splice sites in the e1 region to be joined with the counted splice sites
# from regtools and add read counts for each junction

count_splices <- function(count_data){
  inner_join(tot_unique_spliced_juncs,
             single_junc_stats,
             by = c("junc_ss", "junc_ts")) %>% 
  select(timepoint, chrom, region, transcript_id, strand.x,
         strand.y, junc_name = junc_name.x, junc_ss, junc_ts, read_count, start, end) %>%
  group_by(timepoint, junc_ss, junc_ts, junc_name, .drop = F) %>% 
  reframe(total_read_count = sum(read_count),
          region = region) %>%
  split(.$timepoint) %>%
  map(mutate, tot_reads_time = sum(total_read_count)) %>% 
  do.call("rbind", .) %>%
  mutate(junc_proportion = (total_read_count / tot_reads_time) * 100) %>%
  group_by(region) %>%
  mutate(region_proportion = (sum(total_read_count) / tot_reads_time) * 100) %>% 
  as_tibble()
}

single_count_splice <- count_splices(single_junc_stats)

bulk_count_splice <- count_splices(bulk_junc_stats)

## visualize
single_count_splice %>%
  filter(junc_name == "junc_16") %>% 
  ggplot(aes(timepoint, total_read_count, fill = region)) +
  geom_col(position = "dodge") +
  facet_wrap(~timepoint, scales = "free") +
  scale_fill_igv()

bulk_count_splice %>%
  ggplot(aes(timepoint, total_read_count, fill = region)) +
  geom_col(position = "dodge") +
  facet_wrap(~timepoint, scales = "free") +
  scale_fill_igv()


