#!/usr/bin/env Rscript

library(magrittr)
library(scales)
library(plyr)
library(rtracklayer)
library(glue)
library(ballgown)
library(devtools)
library(genefilter)
library(tidyverse)

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
  select(-ends_with("_1"), -starts_with("start"), -starts_with("end")) %>%
  mutate(trxpt_id = paste0("TRXPT_", seq(1, 28, 1)))

# extract junction pairs (start and end) to make long format
get_junc_pair <- function(j_ss, j_ts){
  return(data.frame(transcript_id = all_trxptome_ss$transcript_id,
                    trxpt_id = all_trxptome_ss$trxpt_id,
                    strand = all_trxptome_ss$strand,
                    region = all_trxptome_ss$region,
                    junc_ss = all_trxptome_ss %>% select(all_of(j_ss)),
                    junc_ts = all_trxptome_ss %>% select(all_of(j_ts))
                    ) %>%
           set_colnames(c("transcript_name", "trxpt_id", "strand", "region",
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
  select(transcript_name, trxpt_id, strand, region, junc_name, junc_ss, junc_ts)


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
    select(timepoint, region = region.x, transcript_name, trxpt_id, strand, junc_name,
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


# ====================================================================
# TRANSCRIPT ABUNDANCE ANALYSIS WITH BALLGOWN
# ====================================================================

# make experimental data dataframe
exp_info <- tribble(~sample_name, ~timepoint, ~ replicate,
                    "abund_12hrsS1", "12h.p.i", "Rep1",
                    "abund_12hrsS3", "12h.p.i", "Rep3",
                    "abund_24hrsS1", "24h.p.i", "Rep1",
                    "abund_24hrsS2", "24h.p.i", "Rep2",
                    "abund_24hrsS3", "24h.p.i", "Rep3",
                    "abund_4hrsS1", "4h.p.i", "Rep1",
                    "abund_4hrsS2", "4h.p.i", "Rep2",
                    "abund_72hrsS1", "72h.p.i", "Rep1",
                    "abund_72hrsS2", "72h.p.i", "Rep2",
                    "abund_72hrsS3", "72h.p.i", "Rep3") %>% data.frame()

# load in expression level files for making ballgown object
bam_files <- list.files("results/hisat2",
                        pattern = "^thev_subset_\\d+hrsS\\d\\.bam$",
                        full.names = T) %>%
  .[c(1:7, 9:11)]

raw_data <- ballgown(dataDir = "results/ballgown",
                     samplePattern = "abund_",
                     pData = exp_info,
                     bamfiles = bam_files)


# extract transcript expression level data
t_exp_levels <- texpr(raw_data, meas = "all") %>%
  as_tibble() %>%
  mutate(trxpt_id = paste0("TRXPT_", seq(1, 28, 1))) %>% 
  pivot_longer(cols = starts_with("FPKM"),
               names_to = "samples",
               values_to = "fpkm") %>%
  mutate(timepoint = case_when(str_detect(samples, "_4hrs") ~ "4h.p.i",
                               str_detect(samples, "_12hrs") ~ "12h.p.i",
                               str_detect(samples, "_24hrs") ~ "24h.p.i",
                               str_detect(samples, "_72hrs") ~ "72h.p.i",
                               TRUE ~ NA_character_),
         timepoint = factor(timepoint,
                            levels = c("4h.p.i", "12h.p.i", "24h.p.i", "72h.p.i"))) %>%
  select(c("timepoint", "t_name", trxpt_id, region = "gene_name", "strand", "start", "end", "fpkm")) %>%
  split(.$timepoint) %>%
  map(mutate, tot_fpkm_tp = sum(fpkm)) %>% 
  do.call("rbind", .)


# ----------------
# expression levels by region
t_exp_lev_byregion <- t_exp_levels %>%
  group_by(timepoint, region, .drop = F) %>%
  reframe(fpkm_reg = sum(fpkm),
          tot_fpkm_tp = mean(tot_fpkm_tp),
          fpkm_reg_percent = round((fpkm_reg / tot_fpkm_tp) * 100, 2)) %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("4h.p.i", "12h.p.i", "24h.p.i", "72h.p.i")))

## ---------------------
# expression levels by transcript
t_expr_each <- t_exp_levels %>%
  mutate(trxpt_id = glue("{region}: {trxpt_id}")) %>%
  group_by(timepoint, t_name, trxpt_id) %>%
  reframe(fpkm_trxpt = sum(fpkm),
          tot_fpkm_tp = mean(tot_fpkm_tp),
          fpkm_trxpt_percent = round((fpkm_trxpt / tot_fpkm_tp) * 100, 2)) %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("4h.p.i", "12h.p.i", "24h.p.i", "72h.p.i")),
         t_name = dplyr::case_match(t_name,
                                    "trxptA_novel" ~ "E1: Hyd_iso_1",
                                    "trxptB_ORF1" ~ "E1: ORF1_novel_iso",
                                    "trxptC_hyd" ~ "E1: Hyd_iso_2",
                                    "trxptD_ORF4" ~ "E1: ORF4_novel",
                                    "DBP" ~ "E2A: DBP",
                                    "ptp_pol1" ~ "E2B: pTP/Pol_iso_1",
                                    "ptp_pol2" ~ "E2B: pTP/Pol_iso_2",
                                    "trunc_hypothetical" ~ "E2B: hypo_trunc",
                                    "IVa2" ~ "IM: IVa2",
                                    "100K:Fib" ~ "E3: 100K->Fiber",
                                    "33K_pVIII_E3_Fib" ~ "E3: 33K->Fiber",
                                    "33K_pVIII_E3_lngr_tts" ~ "E3: 33K->E3_iso_1",
                                    "33K_pVIII_E3_xtra_exon" ~ "E3: 33K->E3_iso_2",
                                    "22K:Fib" ~ "E3: 22K->Fiber",
                                    "E3_Fib_ORF7" ~ "E3: E3->ORF7",
                                    "ORF8" ~ "E4: ORF8",
                                    "33K:E3" ~ "MLP: 33K->E3",
                                    "52K" ~ "MLP: 52K",
                                    "Hexon" ~ "MLP: Hexon",
                                    "Hexon:protease" ~ "MLP: Hexon/Protease",
                                    "TPL_exons_trunc" ~ "MLP: truncated_TPL",
                                    "MLP_Fib_ORF7" ~ "MLP: Fiber/ORF7",
                                    "pIIIa:pVII" ~ "MLP: pIIIa->pVII",
                                    "pVII:protease" ~ "MLP: pVII->Protease_iso_1",
                                    "pVII:protease_2" ~ "MLP: pVII->Protease_iso_2",
                                    "III_pVII" ~ "MLP: III/pVII",
                                    "pVII" ~ "MLP: pVII",
                                    "pX:Hexon" ~ "MLP: pX->Hexon",
                                    .default =  t_name
         )) %>%
  mutate(t_name = factor(t_name,
                         levels = c("E1: Hyd_iso_1", "E1: ORF1_novel_iso",
                                    "E1: Hyd_iso_2", "E1: ORF4_novel",
                                    "E2A: DBP",  "E2B: pTP/Pol_iso_1",
                                    "E2B: pTP/Pol_iso_2", "E2B: hypo_trunc",
                                    "IM: IVa2", "E3: 100K->Fiber",
                                    "E3: 33K->Fiber", "E3: 33K->E3_iso_1",
                                    "E3: 33K->E3_iso_2", "E3: 22K->Fiber",
                                    "E3: E3->ORF7", "E4: ORF8", "MLP: 33K->E3",
                                    "MLP: 52K", "MLP: Hexon", "MLP: Hexon/Protease",
                                    "MLP: truncated_TPL", "MLP: Fiber/ORF7",
                                    "MLP: pIIIa->pVII", "MLP: pVII->Protease_iso_1",
                                    "MLP: pVII->Protease_iso_2", "MLP: III/pVII",
                                    "MLP: pVII", "MLP: pX->Hexon")))


# =================================================================
# REGION-BY-REGION ANALYSES TABLES/DATA
# =================================================================

reg_brkdown_tab <- function(reg){
  
  cleanup_cols <- function(df){
    # cleanup column of duplicates and set them to vectors
    for(col in seq_along(df)){
      
      df[[col]] <- map(df[[col]], ~unique(.x))
      df[[col]] <- map_chr(df[[col]], ~base::paste(.x, collapse = ", "))
    }
    return (df)
  }
  
  reg_ss <- rbind(get_junc_pair("j1_ss", "j1_ts"), get_junc_pair("j2_ss", "j2_ts"),
              get_junc_pair("j3_ss", "j3_ts"), get_junc_pair("j4_ss", "j4_ts"),
              get_junc_pair("j5_ss", "j5_ts"), get_junc_pair("j6_ss", "j6_ts")) %>%
    drop_na() %>%
    mutate(region = factor(region, levels = c("IM", "E1", "E2", "E3", "E4", "MLP")),
           transcript_name = dplyr::case_match(transcript_name,
                                      "trxptA_novel" ~ "Hyd_iso_1",
                                      "trxptB_ORF1" ~ "ORF1_novel_iso",
                                      "trxptC_hyd" ~ "Hyd_iso_2",
                                      "trxptD_ORF4" ~ "ORF4_novel",
                                      "DBP" ~ "DBP",
                                      "ptp_pol1" ~ "pTP/Pol_iso_1",
                                      "ptp_pol2" ~ "pTP/Pol_iso_2",
                                      "trunc_hypothetical" ~ "novel_truncated",
                                      "IVa2" ~ "IVa2",
                                      "100K:Fib" ~ "100K->Fiber",
                                      "33K_pVIII_E3_Fib" ~ "33K->Fiber",
                                      "33K_pVIII_E3_lngr_tts" ~ "33K->E3_iso_1",
                                      "33K_pVIII_E3_xtra_exon" ~ "33K->E3_iso_2",
                                      "22K:Fib" ~ "22K->Fiber",
                                      "E3_Fib_ORF7" ~ "E3->ORF7",
                                      "ORF8" ~ "ORF8",
                                      "33K:E3" ~ "33K->E3",
                                      "52K" ~ "52K",
                                      "Hexon" ~ "Hexon",
                                      "Hexon:protease" ~ "Hexon/Protease",
                                      "TPL_exons_trunc" ~ "TPL",
                                      "MLP_Fib_ORF7" ~ "Fiber/ORF7",
                                      "pIIIa:pVII" ~ "pIIIa->pVII",
                                      "pVII:protease" ~ "pVII->Protease_iso_1",
                                      "pVII:protease_2" ~ "pVII->Protease_iso_2",
                                      "III_pVII" ~ "III/pVII",
                                      "pVII" ~ "pVII",
                                      "pX:Hexon" ~ "pX->Hexon",
                                      .default =  transcript_name)
           ) %>%
    arrange(region) %>%
    group_by(junc_ss, junc_ts) %>%
    reframe(t_name = list(transcript_name),
            trxpt_id = list(trxpt_id),
            region = list(region),
            strand = list(strand)) %>%
    select(t_name, trxpt_id, region, strand, junc_ss, junc_ts)
  
  # unlist columns
  cln_reg_ss <- cleanup_cols(reg_ss)
  
  
  # add junction count data
  cln_reg_ss_counts <- cln_reg_ss %>%
    mutate(junc_ss = as.numeric(junc_ss),
           junc_ts = as.numeric(junc_ts),
           intron_len = (junc_ts - junc_ss) + 1) %>% 
    left_join(., unq_bulk_juncs, by = join_by(junc_ss == start, junc_ts == end)) %>%
    select(-c(trxpts, region.y, accum_tally_j, exact_ss)) %>%
    group_by(timepoint, junc_ss, junc_ts) %>%
    reframe(sum_reads = sum(tot_rds_j_time),
            trxpt_id = list(trxpt_id),
            region = list(region.x),
            strand = list(strand),
            splice_site = list(splice_site),
            intron_len = list(intron_len),
            coding_potential = list(t_name)
            )
  
  cln_reg_ss_counts <- cleanup_cols(cln_reg_ss_counts) %>%
    pivot_wider(names_from = timepoint, values_from = sum_reads) %>%
    dplyr::filter(str_detect(region, reg)) %>%
    dplyr::select(trxpt_id, start = junc_ss, end = junc_ts, splice_site, intron_len,
                  region, strand, "4hpi", "12hpi", "24hpi", "72hpi", coding_potential) %>%
    mutate("4hpi" = ifelse(is.na(`4hpi`), 0, `4hpi`),
           "12hpi" = ifelse(is.na(`12hpi`), 0, `12hpi`),
           "24hpi" = ifelse(is.na(`24hpi`), 0, `24hpi`),
           "72hpi" = ifelse(is.na(`72hpi`), 0, `72hpi`),
           intron_len = glue("{intron_len}bp"))
  
  return(cln_reg_ss_counts)
}
  

# E1 region table analysis for discussion
reg_e1_brkdown <- reg_brkdown_tab("E1")
  
# E2 region table analysis for discussion
reg_e2_brkdown <- reg_brkdown_tab("E2")

# E3 region table analysis for discussion
reg_e3_brkdown <- reg_brkdown_tab("E3")

# E4 region table analysis for discussion
reg_e4_brkdown <- reg_brkdown_tab("E4")


# IM region table analysis for discussion
reg_im_brkdown <- reg_brkdown_tab("IM")

# MLP region table analysis for discussion
reg_mlp_brkdown <- reg_brkdown_tab("MLP")



