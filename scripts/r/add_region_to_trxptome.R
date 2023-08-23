#!/usr/bin/env Rscript

library(tidyverse)


 updated_trxptome_gtf <- read_tsv("results/gffcompare/gffcomp_alltimes.combined.gtf",
              col_names = F, show_col_types = F) %>%
  separate_wider_delim(X9, delim = "; ", names = c("t_id", "others"),
                       too_many = "merge") %>%
  mutate(region = case_when(t_id %in% paste0("transcript_id \"TCONS_000000", c('15', '16', '04', '05'), "\"") ~ "gene_name \"E1\"",
                            t_id %in% paste0("transcript_id \"TCONS_000000", c('03', '17', '18', '28'), "\"") ~ "gene_name \"E2\"",
                            t_id %in% paste0("transcript_id \"TCONS_000000", c('26', '27', '12', '10', '11', '02'), "\"") ~ "gene_name \"E3\"",
                            t_id == "transcript_id \"TCONS_00000014\"" ~ "gene_name \"E4\"",
                            t_id == "transcript_id \"TCONS_00000013\"" ~ "gene_name \"IM\"",
                            TRUE ~ "gene_name \"MLP\""),
         
         t_id = case_when(t_id == "transcript_id \"TCONS_00000015\"" ~ "transcript_id \"trxptB_ORF1\"",
                          t_id == "transcript_id \"TCONS_00000016\"" ~ "transcript_id \"trxptC_hyd\"",
                          t_id == "transcript_id \"TCONS_00000004\"" ~ "transcript_id \"trxptA_novel\"",
                          t_id == "transcript_id \"TCONS_00000005\"" ~ "transcript_id \"trxptD_ORF4\"",
                          # E2
                          t_id == "transcript_id \"TCONS_00000003\"" ~ "transcript_id \"DBP\"",
                          t_id == "transcript_id \"TCONS_00000017\"" ~ "transcript_id \"ptp_pol2\"",
                          t_id == "transcript_id \"TCONS_00000018\"" ~ "transcript_id \"ptp_pol1\"",
                          t_id == "transcript_id \"TCONS_00000028\"" ~ "transcript_id \"trunc_hypothetical\"",
                          # E3
                          t_id == "transcript_id \"TCONS_00000026\"" ~ "transcript_id \"33K_pVIII_E3_Fib\"",
                          t_id == "transcript_id \"TCONS_00000027\"" ~ "transcript_id \"E3_Fib_ORF7\"",
                          t_id == "transcript_id \"TCONS_00000012\"" ~ "transcript_id \"100K:Fib\"",
                          t_id == "transcript_id \"TCONS_00000010\"" ~ "transcript_id \"33K_pVIII_E3_xtra_exon\"",
                          t_id == "transcript_id \"TCONS_00000011\"" ~ "transcript_id \"33K_pVIII_E3_lngr_tts\"",
                          t_id == "transcript_id \"TCONS_00000002\"" ~ "transcript_id \"22K:Fib\"",
                          # E4
                          t_id == "transcript_id \"TCONS_00000014\"" ~ "transcript_id \"ORF8\"",
                          # IM
                          t_id == "transcript_id \"TCONS_00000013\"" ~ "transcript_id \"IVa2\"",
                          # MLP
                          t_id == "transcript_id \"TCONS_00000022\"" ~ "transcript_id \"MLP_Fib_ORF7\"",
                          t_id == "transcript_id \"TCONS_00000020\"" ~ "transcript_id \"pVII:protease\"",
                          t_id == "transcript_id \"TCONS_00000021\"" ~ "transcript_id \"Hexon:protease\"",
                          t_id == "transcript_id \"TCONS_00000019\"" ~ "transcript_id \"pVII:protease_2\"",
                          t_id == "transcript_id \"TCONS_00000006\"" ~ "transcript_id \"52K\"",
                          t_id == "transcript_id \"TCONS_00000007\"" ~ "transcript_id \"33K:E3\"",
                          t_id == "transcript_id \"TCONS_00000008\"" ~ "transcript_id \"Hexon\"",
                          t_id == "transcript_id \"TCONS_00000025\"" ~ "transcript_id \"pVII\"",
                          t_id == "transcript_id \"TCONS_00000023\"" ~ "transcript_id \"pIIIa:pVII\"",
                          t_id == "transcript_id \"TCONS_00000024\"" ~ "transcript_id \"III_pVII\"",
                          t_id == "transcript_id \"TCONS_00000001\"" ~ "transcript_id \"TPL_exons_trunc\"",
                          t_id == "transcript_id \"TCONS_00000009\"" ~ "transcript_id \"pX:Hexon\"",
                          TRUE ~ t_id)
         ) %>% 
  select(X1:X8, region, t_id, others) %>%
  unite("X9", region:others, sep = "; ")

 
 write.table(updated_trxptome_gtf, "results/gffcompare/updated_alltimes.combined.gtf",
             col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
 