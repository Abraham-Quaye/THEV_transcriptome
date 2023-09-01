#!/usr/bin/env Rscript

library(magrittr)
library(glue)
library(tidyverse)
library(scales)

# -------------
# single_jun_files <- list.files("results/hisat2", full.names = T,
#                                pattern = "\\w+\\.bed") %>%
#   setNames(c("t12s1","t12s3", "t24s1", "t24s2", "t24s3","t4s1",
#              "t4s2", "t4s3", "t72s1", "t72s2", "t72s3"))
# # 
# single_junc_stats <- map_dfr(single_jun_files, read_tsv,
#                              col_names = F,
#                              show_col_types = F,
#                              .id = "timepoint") %>%
#   set_colnames(c("timepoint", "chrom", "start", "end", "junc_name",
#                  "read_count", "strand", "thickStart", "thickEnd",
#                  "color", "exons", "exon_sizes", "exon_starts")) %>%
#   arrange(start, end) %>%
#   distinct(start, end, .keep_all = T) %>%
#   # convert from 0-base indexing of .bed format to 1-base indexing
#   mutate(start = start + 1)
# 
# 
# single_sum_stats <- single_junc_stats %>%
#   mutate(timepoint = case_when(timepoint %in% c("t72s1", "t72s2", "t72s3") ~ "72hpi",
#                                timepoint %in% c("t24s1", "t24s2", "t24s3") ~ "24hpi",
#                                timepoint %in% c("t12s1", "t12s3") ~ "12hpi",
#                                timepoint %in% c("t4s1", "t4s2", "t4s3") ~ "4hpi",
#                                .default = timepoint)) %>%
#   group_by(timepoint) %>%
#   summarize(total_junctions = n(),
#             reads10  = sum(read_count > 10),
#             reads100  = sum(read_count > 100),
#             reads1000  = sum(read_count > 1000),
#             total_reads_supporting = sum(read_count),
#             mean_read_supporting = total_reads_supporting/total_junctions,
#             organism = "thev") %>%
#   mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>%
#   arrange(total_junctions)

# ============================BULK=============================================
# =============================================================================
bulk_jun_files <- list.files("results/hisat2/bulk", full.names = T,
                             pattern = "annot_\\d+hrsSS\\.txt") %>% 
  setNames(c("12hpi", "24hpi", "4hpi", "72hpi"))


bulk_junc_stats <- map_dfr(bulk_jun_files, read_tsv,
                           col_types = c("ciicicciiiciiiccc"),
                           show_col_types = F,
                           .id = "timepoint") %>% 
  dplyr::rename(read_count = score) %>%
  arrange(start, end) %>%
  mutate(exact_ss = str_replace(splice_site,
                                "[A-Z]([A-Z])-([A-Z])[A-Z]",
                                "\\1-\\2")) %>%
  mutate(region = case_when(strand == "+" & start > 0 & end <= 2325 ~ "E1",
                            strand == "-" & start >= 2334 & end > 3678 & end <= 18752 ~ "E2",
                            strand == "+" & start >= 18230 & end <= 25168 ~ "E3",
                            strand == "-" & start >= 25191 & end <= 26266 ~ "E4",
                            strand == "-" & start >= 2334 & end <= 3678 ~ "IM",
                            strand == "+" & start < 18230 & end <= 25168 ~ "MLP",
                            TRUE ~ NA)) %>%
  mutate(region = ifelse(is.na(region) & !is.na(gene_names), gene_names, region)
  )

# -----------------------
# function to extract out all unique junctions
extract_unq_juncs <- function(df){
  unqs <- df %>%
  group_by(start, end) %>%
  reframe(Frequency = n(),
          Total_Reads = sum(read_count),
          region = list(region),
          strand = list(strand))

  # cleanup region column of duplicated and NA values
  unqs$region <- map(unqs$region, ~unique(na.omit(.x)))

  # convert region column back to string vector
  unqs$region <- map_chr(unqs$region, ~paste(.x, collapse = ","))

  # remove duplicates from strand column
  unqs$strand <- map(unqs$strand, ~unique(.x))

  # remove ? from strand column
  unqs$strand <- map(unqs$strand, ~.x[.x != "?"])
  
  # convert strand column back to string vector
  unqs$strand <- map_chr(unqs$strand, ~paste(.x, collapse = ","))
  
  unqs <- unqs %>% dplyr::rename(Region = region, Strand = strand)
  return(unqs)
}

total_unq_jncs_submit <- extract_unq_juncs(bulk_junc_stats)
  
# --------------------
# function to extract unique junctions for each timepoint
tot_unique_junc_tp <- function(tp){
  bulk_junc_stats %>%
    filter(timepoint == tp) %>%
    distinct(start, end) %>%
    nrow()
}

tp_unique_juncs <- tibble(total_unq_juncs = c(tot_unique_junc_tp("4hpi"),
                                              tot_unique_junc_tp("12hpi"),
                                              tot_unique_junc_tp("24hpi"),
                                              tot_unique_junc_tp("72hpi")),
                             timepoint = c("4hpi", "12hpi", "24hpi", "72hpi"))

# ------------------------
# calculate overview summary stats
bulk_sum_stats <- bulk_junc_stats %>%
  group_by(timepoint) %>%
  reframe(reads10  = sum(read_count > 10),
          reads100  = sum(read_count > 100),
          reads1000  = sum(read_count > 1000),
          total_reads_tp = sum(read_count),
          organism = "thev") %>% 
  mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>%
  inner_join(., tp_unique_juncs, by = join_by(timepoint)) %>%
  mutate(mean_junc_reads = round(total_reads_tp/total_unq_juncs, 1))


thev_totals <- tibble(sum_junc = sum(bulk_sum_stats$total_unq_juncs),
                      sum_reads10 = sum(bulk_sum_stats$reads10),
                      sum_reads100 = sum(bulk_sum_stats$reads100),
                      sum_reads1000 = sum(bulk_sum_stats$reads1000),
                      mean_mean_read_spt = mean(bulk_sum_stats$mean_junc_reads),
                      sum_junc_reads = sum(bulk_sum_stats$total_reads_tp)) 

# this operation is done here because otherwise, thev_totals will not work
bulk_sum_stats <- bulk_sum_stats %>% 
  mutate(across(where(is.numeric),
                ~ifelse(.x >= 1000000, scientific(.x),
                        ifelse(.x > 1000 & .x < 1000000, comma(.x),
                               .x)
                        )
                )
         ) %>%
  as_tibble()

# ----------------------
# extract how many reads in each bam file
total_reads <- read_table("results/hisat2/coverage/bulk_counts.txt",
                          col_names = T, show_col_types = F)

# ---------------------------
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
  dplyr::rename("organism" = "#rname") %>% 
  mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>% 
  map_at(c(3:10), as.numeric) %>% 
  as_tibble() %>%
  mutate(organism = ifelse(organism == "AY849321.1", "thev", "m.gallopavo")) %>% 
  group_by(organism, timepoint) %>%
  reframe(total_mapped = sum(numreads),
          mean_depth = mean(meandepth), mean_cov = mean(coverage)) %>% 
  cbind(total_reads) %>% 
  mutate(perc_mapped = round((total_mapped / total_reads) * 100, 4)) %>%
  map_at(c(3:7), format, big.mark = ",", digits = 3, drop0trailing = T) %>%
  as_tibble()

join_stats <- left_join(all_covs, bulk_sum_stats, by = c("timepoint", "organism")) %>%
  select(organism, timepoint, total_reads, total_mapped, perc_mapped, mean_depth,
         mean_cov, mean_junc_reads, total_unq_juncs, total_reads_tp, reads10, reads100,
         reads1000) %>%
  as_tibble()

# ------------------
# For Table 1
tab1 <- tibble(Metric = c("Total reads",
                          "Mapped \n (Host)",
                          "Mapped \n (THEV)",
                          "Mean Per Base \n Coverage/Depth",
                          "Total unique \n splice junctions",
                          "Junction coverage \n Total (at least 1 read)",
                          "Junction coverage \n Mean reads",
                          "Junction coverage \n (at least 10 reads)",
                          "Junction coverage \n (at least 100 reads)",
                          "Junction coverage \n (at least 1000 reads)"),
              "4h.p.i" = c(all_covs$total_reads[1],
                           glue("{all_covs$total_mapped[1]} \n ({all_covs$perc_mapped[1]}%)"),
                           glue("{all_covs$total_mapped[5]} \n ({all_covs$perc_mapped[5]}%)"),
                           as.numeric(join_stats$mean_depth[5]),
                           join_stats$total_unq_juncs[5],
                           join_stats$total_reads_tp[5],
                           join_stats$mean_junc_reads[5],
                           join_stats$reads10[5],
                           join_stats$reads100[5],
                           join_stats$reads1000[5]),
              "12h.p.i" = c(all_covs$total_reads[2],
                            glue("{all_covs$total_mapped[2]} \n ({all_covs$perc_mapped[2]}%)"),
                            glue("{all_covs$total_mapped[6]} \n ({all_covs$perc_mapped[6]}%)"),
                            as.numeric(join_stats$mean_depth[6]),
                            join_stats$total_unq_juncs[6],
                            join_stats$total_reads_tp[6],
                            join_stats$mean_junc_reads[6],
                            join_stats$reads10[6],
                            join_stats$reads100[6],
                            join_stats$reads1000[6]),
              "24h.p.i" = c(all_covs$total_reads[3],
                            glue("{all_covs$total_mapped[3]} \n ({all_covs$perc_mapped[3]}%)"),
                            glue("{all_covs$total_mapped[7]} \n ({all_covs$perc_mapped[7]}%)"),
                            join_stats$mean_depth[7],
                            join_stats$total_unq_juncs[7],
                            join_stats$total_reads_tp[7],
                            join_stats$mean_junc_reads[7],
                            join_stats$reads10[7],
                            join_stats$reads100[7],
                            join_stats$reads1000[7]),
              "72h.p.i" = c(all_covs$total_reads[4],
                            glue("{all_covs$total_mapped[4]} \n ({all_covs$perc_mapped[4]}%)"),
                            glue("{all_covs$total_mapped[8]} \n ({all_covs$perc_mapped[8]}%)"),
                            join_stats$mean_depth[8],
                            join_stats$total_unq_juncs[8],
                            join_stats$total_reads_tp[8],
                            join_stats$mean_junc_reads[8],
                            join_stats$reads10[8],
                            join_stats$reads100[8],
                            join_stats$reads1000[8]),
              Total = c(sum(as.numeric(all_covs$total_reads[1:4])),
                        sum(as.numeric(all_covs$total_mapped[1:4])),
                        sum(as.numeric(all_covs$total_mapped[5:8])),
                        sum(as.numeric(str_replace(join_stats$mean_depth[5:8], ",", ""))),
                        total_unq_jncs_submit %>% nrow(),
                        thev_totals$sum_junc_reads,
                        thev_totals$mean_mean_read_spt,
                        thev_totals$sum_reads10,
                        thev_totals$sum_reads100,
                        thev_totals$sum_reads1000)) %>% 
  mutate(Total = ifelse(Total >= 1000000, scientific(Total),
                        ifelse(Total > 1000 & Total < 1000000, comma(Total),
                               Total)))


# ---------------------
# For Table 2A-C
find_sig_juncs <- function(tp){
  total_tp_juncs <- bulk_junc_stats %>% filter(timepoint == tp) %>% pull(read_count) %>% sum()
  
  meta_dt <-  bulk_junc_stats %>%
    filter(timepoint == tp) %>%
    select(timepoint, strand, start, end, splice_site,
           exact_ss, region) %>%
    distinct(start, end, .keep_all = T)
  
  est <- bulk_junc_stats %>%
  filter(timepoint == tp) %>%
  group_by(timepoint, start, end) %>%
  reframe(tot_tp_j_rds = sum(read_count)
          ) %>%
  mutate(perc = round((tot_tp_j_rds/ total_tp_juncs) * 100, 1),
         intron_len = end - start) 
  
  return(
    left_join(meta_dt, est, by = join_by(timepoint, start, end)) %>%
    filter(perc >= 1) %>%
    arrange(desc(perc)) %>%
    distinct(start, end, .keep_all = T) %>%
    mutate(tot_tp_j_rds = comma(tot_tp_j_rds),
           intron_len = comma(intron_len),
           counts_percent = glue("{tot_tp_j_rds} ({perc}%)")) %>%
    select(-c(perc)) %>%
    mutate(intron_len = glue("{intron_len} bp"),
           strand = ifelse(strand == "+", "\u002B", "\u002D")) %>% 
    dplyr::rename(Start = start, End = end, Timepoint = timepoint,
                  Strand = strand, Region = region, Splice_Site = splice_site,
                  "Splice \n Acceptor-Donor" = exact_ss, Reads = tot_tp_j_rds,
                  "Intron Length" = intron_len, Reads_Percentage = counts_percent)
  )
}
    
sig_12_juncs <- find_sig_juncs("12hpi")

sig_24_juncs <- find_sig_juncs("24hpi")

sig_72_juncs <- find_sig_juncs("72hpi")

# ----------------------
# For Supplementary Tables 1A-C

tp_reg_expr_abund <- function(tp){
  total_tp_rds <- bulk_junc_stats %>% filter(timepoint == tp) %>% pull(read_count) %>% sum()
  
  tp_reg_stats <- bulk_junc_stats %>%
    filter(timepoint == tp) %>%
    group_by(timepoint, region) %>%
    reframe(strand = list(strand),
            junc_freq = n(),
            sum_reads_reg = sum(read_count),
            perc_reg = glue("{round((sum_reads_reg/ total_tp_rds) * 100, 1)}%")) %>%
    dplyr::arrange(desc(sum_reads_reg))
  
  # remove duplicates from strand column
  tp_reg_stats$strand <- map(tp_reg_stats$strand, ~unique(.x))
  
  # remove ? from strand column
  tp_reg_stats$strand <- map(tp_reg_stats$strand, ~.x[.x != "?"])
  
  # convert strand column back to string vector
  tp_reg_stats$strand <- map_chr(tp_reg_stats$strand, ~paste(.x, collapse = ","))
  
  return(tp_reg_stats)
}

tp_reg_expr_12 <- tp_reg_expr_abund("12hpi")

tp_reg_expr_24 <- tp_reg_expr_abund("24hpi")

tp_reg_expr_72 <- tp_reg_expr_abund("72hpi")

