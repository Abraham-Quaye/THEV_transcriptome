#!/usr/bin/env Rscript

library(magrittr)
library(scales)
library(tidyverse)
library(plyr)
library(rtracklayer)
library(glue)
library(ggsci)
library(patchwork)

# load in data from from full transcriptome to access the splice sites
all_trxptome_ss <- import("results/gffcompare/updated_alltimes.combined.gtf") %>%
  as_tibble() %>%
  select(chr = seqnames, gene_id, transcript_id, gene_name, strand, type,
         start, end, orf_size = width, exon_number) %>%
  filter(type != "transcript") %>%
  mutate(type = glue("{type}_{exon_number}")) %>%
  group_by(transcript_id, strand, region = gene_name) %>%
  reframe(start = list(start),
          end = list(end)) %>%
  unnest_wider(c(start, end), names_sep = "_") %>%
  arrange(region) %>%
  mutate(j1_ss = end_1, j1_ts = start_2,
         j2_ss = end_2, j2_ts = start_3,
         j3_ss = end_3, j3_ts = start_4,
         j4_ss = end_4, j4_ts = start_5,
         j5_ss = end_5, j5_ts = start_6,
         j6_ss = end_6, j6_ts = start_7
  ) %>% 
  select(-ends_with("_1"), -starts_with("start"), -starts_with("end"))


# extract junction pairs (start and end) to make long format
get_junc_pair <- function(j_ss, j_ts){
  return(data.frame(transcript_id = all_trxptome_ss$transcript_id,
             strand = all_trxptome_ss$strand,
             region = all_trxptome_ss$region,
             junc_ss = all_trxptome_ss %>% select(all_of(j_ss)),
             junc_ts = all_trxptome_ss %>% select(all_of(j_ts))
             ) %>%
           set_colnames(c("transcript_id", "strand", "region",
                                  "junc_ss", "junc_ts")) %>%
           as_tibble()
  )
}

# Unique splice junctions in the finalized transcriptome
unq_trxptome_juncs <- rbind(get_junc_pair("j1_ss", "j1_ts"), get_junc_pair("j2_ss", "j2_ts"),
                                  get_junc_pair("j3_ss", "j3_ts"), get_junc_pair("j4_ss", "j4_ts"),
                                  get_junc_pair("j5_ss", "j5_ts"), get_junc_pair("j6_ss", "j6_ts")) %>%
  drop_na() %>%
  mutate(region = factor(region, levels = c("IM", "E1", "E2", "E3", "E4", "MLP"))) %>%
  arrange(region) %>%
  distinct(junc_ss, junc_ts, .keep_all = T) %>%
  mutate(junc_name = glue("junc_{c(1:nrow(.))}"),
         region = factor(region, levels = c("E1", "E2", "E3", "E4", "IM", "MLP"))) %>%
  select(transcript_id, strand, region, junc_name, junc_ss, junc_ts)


# load in counted splice sites from regtools output to quantify 
# the splices sites of any specified region 
# -----------------

bulk_jun_files <- list.files("results/hisat2/bulk", full.names = T,
                             pattern = "annot_\\d+hrsSS\\.txt") %>% 
  setNames(c("12hpi", "24hpi", "4hpi", "72hpi"))


bulk_junc_stats <- map_dfr(bulk_jun_files, read_tsv,
                           col_types = c("ciicicciiiciiiccc"),
                           show_col_types = F,
                           .id = "timepoint")

assign_region <- function(df){
  df %>%
    dplyr::rename(read_count = score) %>%
    dplyr::arrange(start, end) %>%
    dplyr::mutate(exact_ss = str_replace(splice_site,
                                  "[A-Z]([A-Z])-([A-Z])[A-Z]",
                                  "\\1-\\2")) %>%
    mutate(region = case_when(strand == "+" & start > 0 & end <= 2325 ~ "E1",
                              strand == "-" & start >= 2334 & end > 3678 & end <= 18752 ~ "E2",
                              strand == "+" & start >= 18230 & end <= 25168 ~ "E3",
                              strand == "-" & start >= 25191 & end <= 26266 ~ "E4",
                              strand == "-" & start >= 2334 & end <= 3678 ~ "IM",
                              strand == "+" & start < 18230 & end <= 25168 ~ "MLP",
                              TRUE ~ NA)) %>% 
    mutate(region = ifelse(is.na(region) & !is.na(gene_names), gene_names, region))
}

bulk_junc_stats <- assign_region(bulk_junc_stats)

# na_reg <- bulk_junc_stats %>% filter(is.na(region)) %>% distinct(start, end) %>% nrow()

group_unq_juncs <- function(df){
  df %>%
  group_by(timepoint, start, end, exact_ss, splice_site) %>%
  reframe(accum_tally_j = n(),
          region = list(region),
          tot_rds_j_time = sum(read_count),
          trxpts = list(transcripts))
}

unq_bulk_juncs <- group_unq_juncs(bulk_junc_stats)

cleanup_reg_name <- function(df){
  # cleanup region column of duplicated and NA values
  df$region <- map(df$region, ~unique(na.omit(.x)))
  
  # convert region column back to string vector
  df$region <- map_chr(df$region, ~paste(.x, collapse = ","))
  
  df <- mutate(df, region = ifelse(region == "", NA_character_, region))
  return (df)
}

unq_bulk_juncs <- cleanup_reg_name(unq_bulk_juncs)

# -----------------------
# filter out splice sites in the e1 region to be joined with the counted splice sites
# from regtools and add read counts for each junction

count_trxptome_ss <- function(count_data){
  spl_sites <- left_join(unq_trxptome_juncs,
                         count_data,
                         by = join_by(junc_ss == start, junc_ts == end)) %>%
    select(timepoint, region = region.x, transcript_id, strand, junc_name,
           junc_ss, junc_ts, exact_ss, splice_site, trxpts) %>%
    mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>%
    distinct(junc_name, timepoint, region, .keep_all = T)
  
  
  counts <- left_join(unq_trxptome_juncs,
                      count_data,
                      by = join_by(junc_ss == start, junc_ts == end)) %>%
  select(timepoint, region = region.x, junc_name, accum_tally_j, tot_rds_j_time) %>%
  mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>%
  group_by(timepoint, junc_name, region) %>% 
  reframe(accumul_juncs_reg = sum(accum_tally_j),
          total_rd_j_tp_reg = sum(tot_rds_j_time)) %>%
  split(.$timepoint) %>%
  map(mutate, tot_reads_time = sum(total_rd_j_tp_reg)) %>% 
  do.call("rbind", .) %>%
  mutate(region = factor(region, levels = c("E1", "E2", "E3", "E4", "IM", "MLP"))) %>%
  group_by(timepoint, region, .drop = F) %>% 
  reframe(junc_name = junc_name,
          total_rd_j_tp_reg = total_rd_j_tp_reg,
          tot_reads_time = tot_reads_time,
          count_j_tp_rg = sum(accumul_juncs_reg)) %>%
  arrange(timepoint, region)
  
  return(
    left_join(spl_sites, counts, by = c("junc_name", "timepoint")) %>%
      dplyr::rename(region = region.x) %>%
      arrange(timepoint, region)
  )
}

bulk_count_trxptome_ss <- count_trxptome_ss(unq_bulk_juncs) 

summarize_trxptome_ss <- function(df){
  df_mod <- df %>%
    group_by(timepoint, region, .drop = FALSE) %>%
    reframe(mean_tot_rds_time = mean(tot_reads_time, na.rm = T),
            tot_rds_tp_reg = sum(total_rd_j_tp_reg, na.rm = T),
            reg_junc_percent = round((tot_rds_tp_reg / mean_tot_rds_time) * 100, 2)) %>%
    replace_na(replace = list(mean_tot_rds_time = 0,
                              reg_junc_percent = 0))
  return(df_mod)
}

bulk_count_trxptome_ss <- summarize_trxptome_ss(bulk_count_trxptome_ss)

## visualize
plot_time_genexpression <- function(count_ss){
  up_lim = max(count_ss$reg_junc_percent) %>% plyr::round_any(., 10, f = ceiling)
  count_ss %>%
    filter(timepoint != "4hpi") %>%
    ggplot(aes(timepoint, reg_junc_percent, group = region, color = region)) +
    geom_line(linewidth = 3) +
    geom_point(size = 5) +
    scale_colour_jco() +
    scale_y_continuous(expand = c(0.01, 0.01),
                       breaks = c(seq(0, up_lim, 10)),
                       labels = scales::label_percent(scale = 1)) +
    scale_x_discrete(expand = c(0.005, 0.01)) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(plot.margin = margin(rep(30, 4)),
          panel.grid.major.y = element_line(linewidth = 0.4, color = "grey",
                                            linetype = "dashed"),
          axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 15)),
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 14, face = "bold"),
          legend.justification = c(0, 0),
          legend.position = c(0.05, 0.9),
          legend.key.width = unit(2, "cm"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, face = "bold", colour = "black"),
          legend.text.align = 0) +
    guides(color = guide_legend(label.position = "top",
                                label.hjust = 0.5,
                                label.vjust = 1,
                                direction = "horizontal",
                                nrow = 1,
                                ))
}


# plot and save abundances of juncs found in final trxptome
trxptome_juncs <- plot_time_genexpression(bulk_count_trxptome_ss) +
  labs(x = element_blank(),
       y = "Relative Abundances of Junctions in Final Transcriptome")

# ggsave(plot = trxptome_juncs,
#        filename = "results/r/figures/temporal_gene_expression.png",
#        dpi = 500, width = 14, height = 12)


# prepare data to plot abundances of all juncs extracted from BAM files
plot_all_junc_abunds <- bulk_junc_stats %>%
  group_by(timepoint, region) %>%
  reframe(accum_tally_j = n(),
          total_rd_j_tp_reg = sum(read_count)) %>%
  drop_na(region) %>% 
  split(.$timepoint) %>%
  map(mutate, tot_reads_time = sum(total_rd_j_tp_reg)) %>% 
  do.call("rbind", .) %>%
  mutate(region = factor(region, levels = c("E1", "E2", "E3", "E4", "IM", "MLP"))) %>%
  group_by(timepoint, region, .drop = F) %>% 
  reframe(total_rd_j_tp_reg = total_rd_j_tp_reg,
          tot_reads_time = tot_reads_time,
          count_j_tp_rg = sum(accum_tally_j)) %>%
  arrange(timepoint, region) %>%
  group_by(timepoint, region, .drop = FALSE) %>%
  reframe(mean_tot_rds_time = mean(tot_reads_time, na.rm = T),
          tot_rds_tp_reg = sum(total_rd_j_tp_reg, na.rm = T),
          reg_junc_percent = round((tot_rds_tp_reg / mean_tot_rds_time) * 100, 2)) %>%
  replace_na(replace = list(mean_tot_rds_time = 0,
                            reg_junc_percent = 0)) %>%
  mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi")))

all_juncs <- plot_time_genexpression(plot_all_junc_abunds) +
  labs(x = element_blank(),
       y = "Relative Abundances of All Junctions")

patch_expr <- (all_juncs | trxptome_juncs)

ggsave("junc_abundances.png",
       plot = patch_expr, path = "results/r/figures",
       width = 20, height = 12, dpi = 500)


# -------------------
## splice donor and acceptor frequencies

# acceptors and donors in trxptome
bulk_trxptome_ss_seq <- count_trxptome_ss(unq_bulk_juncs) %$%
  count(exact_ss) %>%
  mutate(sum_junc_count = sum(freq),
         percent_abund = (freq / sum_junc_count) * 100) %>%
  dplyr::rename(ss_seq = x)

# plot trxptome acceptors and donor frequencies
bulk_trxptome_ss_seq %>%
  ggplot(aes(ss_seq, percent_abund, color = ss_seq)) +
  geom_point(show.legend = F, size = 10) +
  geom_segment(aes(x = ss_seq, xend = ss_seq,
                   y = 0, yend = percent_abund),
               linewidth = 2,
               show.legend = F) +
  geom_text(aes(label = glue("{round(percent_abund, 1)}%")),
             nudge_y = 3.5,
            size = 10, fontface = "bold", color = "#000000") +
  labs(x = "Splice Site Donor-Acceptor",
       y = "Frequency") +
  scale_color_lancet() +
  scale_y_continuous(expand = c(0.01,0.01),
                     labels = scales::label_percent(scale = 1)) +
  scale_x_discrete(expand = c(0.13, 0.13)) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(plot.margin = margin(rep(30, 4)),
        panel.grid.major.y = element_line(linewidth = 0.4, color = "grey",
                                          linetype = "dashed"),
        axis.title = element_text(size = 16, face = "bold", margin = margin(r = 15, t = 25)),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(color = "grey40"),
        axis.ticks.length = unit(0, "pt")
        )

# All acceptors and donors
all_ss_seq <- bulk_junc_stats %>%
  mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>% 
  group_by(timepoint, exact_ss) %>%
  reframe(tally = n(),
          total_reads = sum(read_count)) %>%
  split(.$timepoint) %>%
  map(mutate, tot_tally_tp = sum(tally)) %>% 
  do.call("rbind", .) %>%
  mutate(percent_abund = (tally / tot_tally_tp) * 100) %>% 
  ggplot(aes(exact_ss, percent_abund, group = timepoint, color = timepoint)) +
  geom_point(show.legend = F, size = 10) +
  geom_segment(aes(x = exact_ss, xend = exact_ss,
                   y = 0, yend = percent_abund),
               linewidth = 2,
               show.legend = F) +
  geom_text(aes(label = glue("{round(percent_abund, 1)}%")),
            nudge_y = 4.5,
            size = 5, fontface = "bold", color = "#000000") +
  facet_wrap(~ timepoint, scales = "free") +
  labs(x = "Splice Site Donor-Acceptor",
       y = "Frequency") +
  scale_color_lancet() +
  scale_y_continuous(expand = c(0.01,0.01),
                     limits = c(0, 105),
                     labels = scales::label_percent(scale = 1)) +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(plot.margin = margin(rep(30, 4)),
        panel.grid.major.y = element_line(linewidth = 0.4, color = "grey",
                                          linetype = "dashed"),
        axis.title = element_text(size = 16, face = "bold", margin = margin(r = 15, t = 15)),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(color = "grey40"),
        axis.ticks.length = unit(0, "pt"),
        strip.background = element_rect(linewidth = 0.1, fill = "grey"),
        strip.text.x = element_text(size = 12, face = "bold", margin = margin(0.2,0,0.2,0, "cm")),
        strip.placement = "outside",
        panel.spacing = unit(1, "cm")
  )



