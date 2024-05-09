
#!/usr/bin/env Rscript

library(magrittr)
library(glue)
library(rtracklayer)
library(patchwork)
library(ggtext)
library(plyr)
library(tidyverse)
library(ballgown)


# Import GTF and Wrangle for plotting ==========================================
#===============================================================================

spliced_gtf <- import("results/gffcompare/gffcomp_alltimes.combined.gtf") %>%
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
  arrange(strand, region)

trxpt_order <- list(spliced_gtf$transcript_id) %>% unlist()

spliced_gtf <- spliced_gtf %>% 
  arrange(start_1, end_1) %>%
  mutate(trxpt_id = paste0("TRXPT_", seq(1, 28, 1)),
         transcript_id = factor(transcript_id, levels = trxpt_order)) %>%
  arrange(transcript_id) %>%
  # Manually add transcript 29 from initial stringtie results before gffcompare merge
  rbind(tibble(transcript_id = "TCONS_00000029", strand = "+", start_1 = 18230, start_2 = 18230,
               start_3 = 20162, start_4 = NA, start_5 = NA, start_6 = NA, start_7 = NA, start_8 = NA,
               end_1 = 20732, end_2 = 18350, end_3 = 20732, end_4 = NA, end_5 = NA, end_6 = NA,
               end_7 = NA, end_8 = NA, region = "E3", trxpt_id = "TRXPT_29")) %>%
  # Manually add transcript 31 from junction validation discovery
  rbind(tibble(transcript_id = "TCONS_00000031", strand = "-", start_1 = 2334, start_2 = 2334,
               start_3 = 10981, start_4 = 18159, start_5 = 18684, start_6 = NA, start_7 = NA, start_8 = NA,
               end_1 = 18751, end_2 = 7062, end_3 = 11079, end_4 = 18189, end_5 = 18751, end_6 = NA,
               end_7 = NA, end_8 = NA, region = "E2", trxpt_id = "TRXPT_31"))


# trxpt id and names info

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
                    "abund_72hrsS3", "72h.p.i", "Rep3") %>%
  data.frame()

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
  select(trxpt_id, t_name, num_exons) %>%
  # add trxpt_29
  rbind(tibble(trxpt_id = "TRXPT_29", t_name = "22K", num_exons = 2)) %>%
  # add trxpt_31
  rbind(tibble(trxpt_id = "TRXPT_31", t_name = "Adpol", num_exons = 4))


comp_spliced_gtf <- left_join(spliced_gtf, t_exp_levels,
                              by = "trxpt_id") %>%
  mutate(sstcodon = c(# E1
                      rep(211, 3), 1965, 
                      # E3
                      rep(18230, 4), 20769, 22519,
                      # MLP
                      NA, 8570, 20769, 13610, 12712, 12347, 12347, 13610, 22519,
                      9462, 11001, 12347,
                      # E2
                      18186, 10995, 10995, NA, 
                      # E4
                      26246, 
                      # IM
                      3616, 
                      # trxpt_29
                      18230,
                      # trxpt_31
                      6768),
         stpcodon = c(# E1 
                      2312, 1953, rep(2312, 2), 
                      # E3
                      20262, 20699, 20457, 20227, 21371, 23883,
                      # MLP
                      NA, 9472, 21371, 16330, 12888, 12709, 12709, 16330, 23883,
                      10979, 12347, 12709,
                      # E2
                      16973, 6765, 6765, NA, 
                      # E4
                      25204, 
                      # IM
                      2334,
                      # trxpt_29
                      20262,
                      # trxpt_31
                      3430
                      ),
         secSSC = c(rep(NA, 4), 
                    # E3
                    rep(20769, 3), 20142, 21214, 24512, 
                    #MLP
                    NA, NA, 21214, NA, 12906, NA, NA, 16188, 24512, rep(NA, 3),
                    rep(NA, 8)),
         secSTC = c(rep(NA, 4), 
                    # E3
                    rep(21371, 3), 20411, 22116, 25168, 
                    # MLP
                    NA, NA, 22116, NA, 13601, NA, NA, 16976, 25168, rep(NA, 3),
                    rep(NA, 8)))


# load in data and prepare predicted ORFs only

predicted_orfs <- import("results/stringtie/all_merged.gtf") %>% 
  as_tibble() %>% 
  drop_na(reference_id) %>% 
  distinct(start, end, width, .keep_all = T) %>% 
  select(chr = seqnames, start, end, orf_size = width, strand, type, gene_id,
         transcript_id, reference_id, exon_number) %>% 
  group_by(transcript_id, strand) %>%
  reframe(start = list(start),
          end = list(end)) %>%
  unnest_wider(c(start, end), names_sep = "_")

# pull out the orf names
orf_names <- import("results/stringtie/all_merged.gtf") %>% 
  as_tibble() %>%
  drop_na(reference_id) %>% 
  distinct(start, end, width, .keep_all = T) %>% 
  select(gene_name = reference_id, transcript_id) %>%
  distinct(gene_name, .keep_all = T)

# assign orf names back to orf data
predicted_orfs <- inner_join(orf_names, predicted_orfs, by = "transcript_id")

predicted_orfs <- predicted_orfs %>%
  mutate(region = case_when(gene_name %in% c("ORF1", "Hyd") ~ "E1",
                            gene_name %in% c("ORF8") ~ "E4",
                            gene_name %in% c("DBP", "pTP", "AdPOL") ~ "E2",
                            gene_name  %in% c("E3", "100K", "33K", "22K", "pVIII") ~ "E3",
                            gene_name == "IVa2" ~ "IM",
                            gene_name == "UXP" ~ "LONE",
                            TRUE ~ "MLP")) %>%
  arrange(strand, region) %>%
  mutate(start_4 = NA,
         end_4 = NA) %>% 
  # manually added transcripts not from StringTie
  rbind(tibble( gene_name = "ORF4", transcript_id = "gp04", strand = "+",
                start_1 = 1965, start_2 = NA, start_3 = NA, start_4 = NA,
                end_1 = 2312, end_2 = NA, end_3 = NA, end_4 = NA,
                region = "E1")) %>%
  # Full ORF1 transcript
  rbind(tibble( gene_name = "TRXPT_2B", transcript_id = "fORF1", strand = "+",
                start_1 = 54, start_2 = 54, start_3 = 54, start_4 = NA,
                end_1 = 2325, end_2 = 2325, end_3 = 2325, end_4 = NA,
                region = "E1")) %>%
  #TRXPT_21B (DBP isoform2)
  rbind(tibble( gene_name = "TRXPT_21B", transcript_id = "DBP_iso2", strand = "-",
                start_1 = 16934, start_2 = 16934, start_3 = 18684, start_4 = NA,
                end_1 = 18751, end_2 = 18189, end_3 = 18751, end_4 = NA,
                region = "E2")) %>%
  #TRXPT_25B (100K isoform2)
  rbind(tibble( gene_name = "TRXPT_25B", transcript_id = "100K", strand = "+",
                start_1 = 18230, start_2 = NA, start_3 = NA, start_4 = NA,
                end_1 = 23702, end_2 = NA, end_3 = NA, end_4 = NA,
                region = "E3")) %>% 
  #TRXPT_30 (22K long isoform)
  rbind(tibble( gene_name = "TRXPT_30", transcript_id = "l22K", strand = "+",
                start_1 = 18230, start_2 = 18230, start_3 = 18717, start_4 = 20162,
                end_1 = 23884, end_2 = 18350, end_3 = 18768, end_4 = 23884, region = "E3"))

combined_gtf <- plyr::rbind.fill(comp_spliced_gtf, predicted_orfs) %>% 
  arrange(strand, region) %>% as_tibble() %>%
  mutate(trxpt_id = ifelse(is.na(trxpt_id), gene_name, trxpt_id),
         sstcodon = dplyr::case_match(transcript_id,
                                      "fORF1" ~ 211,
                                      "DBP_iso2" ~ 18013,
                                      "100K" ~ 18230,
                                      "l22K" ~ 18230,
                                      .default = sstcodon),
         stpcodon = dplyr::case_match(transcript_id,
                                      "fORF1" ~ 1953,
                                      "DBP_iso2" ~ 16973,
                                      "100K" ~ 20227,
                                      "l22K" ~ 20411,
                                      .default = stpcodon))



plot_full_trxptome <- function(combined_gtf, trxptome_part, trxpt_part2=NULL){
  regs <- c(trxptome_part, trxpt_part2)
  
  # =======================================================
  # TRANSCRIPT X AND Y POSITIONS FOR FULL TRANSCRIPTOME MAP
  # =======================================================
  # add y-axis positions to postive strand E1 transcripts
  ypos_vals <- c(26:45)
  ypos_pos_e1 <- combined_gtf %>% filter(region == "E1")
  ypos_pos_e1$ypos <- rep(ypos_vals, length.out = nrow(ypos_pos_e1))
  
  # add y-axis positions to postive strand MLP transcripts
  ypos_pos_mlp <- combined_gtf %>% filter(region == "MLP")
  ypos_pos_mlp$ypos <- rep(ypos_vals, length.out = nrow(ypos_pos_mlp))
  
  # add y-axis positions to postive strand E3 transcripts
  ypos_pos_e3 <- combined_gtf %>% filter(region == "E3")
  ypos_pos_e3$ypos <- rep(c(26:45), length.out = nrow(ypos_pos_e3))
  
  
  # add y-axis positions to negative strand transcripts
  
  ypos_e4 <- combined_gtf %>% filter(region == "E4")
  ypos_e4$ypos <- rep(c(23, 24), length.out = nrow(ypos_e4))
  
  ypos_e2 <- combined_gtf %>% filter(region == "E2")
  ypos_e2$ypos <- rep(c(24:19), length.out = nrow(ypos_e2))
  
  ypos_im <- combined_gtf %>% filter(region == "IM")
  ypos_im$ypos <- rep(c(21, 21.5), length.out = nrow(ypos_im))
  
  ypos_uxp <- combined_gtf %>% filter(region == "LONE")
  ypos_uxp$ypos <- rep(c(24), length.out = nrow(ypos_uxp))
  
  
  trxpt_regs <- if(trxptome_part == "all"){
    c(unique(combined_gtf$region))  
  } else{c(trxptome_part)}
  
  # recreate full dataframe with y-axis positions
  combined_gtf <- rbind(ypos_pos_e1, ypos_pos_mlp, ypos_pos_e3,
                        ypos_e4, ypos_e2, ypos_im, ypos_uxp) %>% 
    mutate(color = ifelse(!is.na(gene_name), "grey80", "#ff0000"),
           ypos = case_when(gene_name == "Hyd" ~ 30,
                            gene_name == "100K" ~ 25.5,
                            gene_name == "52K" ~ 39,
                            gene_name == "pVII" ~ 36.5,
                            gene_name == "Hexon" ~ 34.5,
                            gene_name == "E3" ~ 26.25,
                            gene_name == "22K" ~ 31.5,
                            gene_name == "pVIII" ~ 31.5,
                            gene_name == "ORF4" ~ 28.5,
                            gene_name == "TRXPT_2B" ~ 29.5,
                            gene_name == "TRXPT_21B" ~ 23.5,
                            gene_name == "pTP" ~ 23.5,
                            gene_name == "DBP" ~ 22.5,
                            gene_name == "AdPOL" ~ 22.5,
                            gene_name == "TRXPT_25B" ~ 28.5,
                            trxpt_id == "TRXPT_25" ~ 29.5,
                            trxpt_id == "TRXPT_26" ~ 30.5,
                            trxpt_id == "TRXPT_29" ~ 32,
                            trxpt_id == "TRXPT_30" ~ 33,
                            trxpt_id == "TRXPT_31" ~ 20.5,
                            TRUE ~ ypos)) %>% 
    mutate(color = ifelse(trxpt_id %in% c("ORF4", "TRXPT_2B",
                                          "TRXPT_21B", "TRXPT_25B",
                                          "TRXPT_30", "TRXPT_31"),
                          "#000000", color)) %>%
    as_tibble() %>%
    filter(region %in% regs)
  
  # ==================================================================
  # TRANSCRIPT POSITIONS FOR REGION-BY-REGION TRANSCRIPT MAPS
  # ==================================================================
  
  # mlp
  if(trxptome_part == "MLP"){
    mlp_adj_ypos <- combined_gtf %>%
      filter(region == trxptome_part & ypos > 38) %>%
      mutate(ypos = rep(c(38, 39), length.out = nrow(.)))
    
    combined_gtf <- combined_gtf %>% filter(!(region == 'MLP' & ypos > 38)) %>%
      rbind(., mlp_adj_ypos)
  }
  
  # ===================================================================
  # SINGLE EXON TRXPTS AND REGION LABELS FOR FULL TRANSCRIPTOME =======
  single_orfs <- combined_gtf %>% filter(is.na(start_2), is.na(end_2))
  
  genome_ruler <- tibble(x = seq(0, 26000, 2000),
                         y = 24.8,
                         yend = 25.3,
                         lab = paste0(seq(0,26, 2), "kb")
  )
  
  
  # ========================================================
  # plot transcripts 
  # ========================================================
  
  # Extract the column names that start with "start" or "end"
  start_cols <- names(combined_gtf)[str_detect(names(combined_gtf), "^start")]
  end_cols <- names(combined_gtf)[str_detect(names(combined_gtf), "^end")]
  
  # 
  x_axis_start <- min(combined_gtf[, start_cols[1]])
  if(!str_detect(trxptome_part, "all|E1")){
    if(trxptome_part == "E4"){
      x_axis_start <- x_axis_start - (x_axis_start * 0.01)
    }else if(trxptome_part == "E2"){
      x_axis_start <- x_axis_start - (x_axis_start * 0.25)
    }else{x_axis_start <- x_axis_start - (x_axis_start * 0.03)}
    x_axis_start <- plyr::round_any(x_axis_start, 100)
  } else{x_axis_start <- 0}
  
  #
  if(trxptome_part != "all"){
    # calc the max for all end colums
    x_axis_candidates <- c()
    for(col in seq_along(end_cols)){
      x_axis_candidates[col] <- max(combined_gtf[, end_cols[col]], na.rm = T)
    }
    # pull the max of the candidate maximums
    x_axis_end <- max(x_axis_candidates, na.rm = T)
    # modify for use in axis
    if(trxptome_part == "E4"){
      x_axis_end <- x_axis_end + (x_axis_end * 0.01)
    }else{x_axis_end <- x_axis_end + (x_axis_end * 0.03)}
    x_axis_end <- plyr::round_any(x_axis_end, 100, f = ceiling)
    x_axis_end <- ifelse(x_axis_end > 26266, 26266, x_axis_end)
  } else{x_axis_end <- 26266}
  
  
  if(trxptome_part == "all"){
    y_lims <- c(NA, 48)
  } else if(trxptome_part == "E4"){
    y_lims <- c(22.5, NA)
  }else if(trxptome_part == "E2"){
    y_lims <- c(20, NA)
  }else if(trxptome_part == "IM"){
    y_lims <- c(20.5, NA)
  }else{y_lims <- c(NA, NA)}
  
  splice_map <- combined_gtf %>%
    ggplot() +
    # line representing whole genome
    geom_segment(aes(x = x_axis_start, xend = x_axis_end, y = 25, yend = 25),
                 linewidth = 2, color = "black") +
    
    # plot full trxpts: start_pos to end_pos
    geom_segment(aes_string(x = start_cols[1], xend = end_cols[1],
                            y = "ypos", yend = "ypos"),
                 linetype = "dotted", color = combined_gtf$color) +
    
    # plot single orfs
    geom_segment(data = single_orfs, aes(x = start_1, xend = end_1, y = ypos, yend = ypos),
                 linewidth = 7, color = single_orfs$color) +
    
    # aesthetics
    scale_x_continuous(expand = c(0.01,0.01),
                       limits = c(x_axis_start, (x_axis_end + (x_axis_end * 0.005))),
                       breaks = seq(0, 26000, 2000), 
                       labels = paste0(seq(0,26, 2), "kb")) +
    scale_y_continuous(expand = c(0.01,0.01),
                       limits = y_lims) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(rep(20, 4)),
          plot.background = element_rect(fill = "#ffffff", color = c("grey")),
          panel.background = element_rect(fill = "#ffffff"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.title = element_blank(),
          axis.text= element_blank(),
          axis.line.x = element_blank(),
          axis.ticks = element_blank(),
    )
  
  # Loop through the remaining exons (start and end columns) and 
  # add a geom_segment layer for each pair
  for (i in seq_along(start_cols)[-1]) {
    if(!base::is.na(start_cols[i])){
      splice_map <- splice_map +
        geom_segment(aes_string(x = start_cols[i], xend = end_cols[i], y = "ypos", yend = "ypos"),
                     color = combined_gtf$color, linewidth = 7)
    }
  }
  
  # Add labels for full plot -----------------
  splice_map <- splice_map +
    # label each transcript
    geom_richtext(aes(x = end_1, y = ypos, label = trxpt_id,
                      hjust = 0, fontface = "bold"), size = 3,
                  label.size = NA, label.padding = unit(0, "pt"),
                  fill = NA, label.margin = margin(l = 5),
                  show.legend = F)
  
  return(splice_map)
}


# ===============================================================
# REGION-BY-REGION BREAKDOWN PLOTS
# ===============================================================

brkdown_reg_plots <- function(reg, reg2=NULL){
  
  # genome ruler for region-by-region breakdown
  if(reg %in% c("E1", "IM")){
    genome_ruler <- tibble(x = seq(0, 26000, 200),
                           y = 24.8,
                           yend = 25.3,
                           lab = paste0(seq(0, 26, 0.2), "kb")
    )
  }else if(reg %in% c("E2", "MLP")){
    genome_ruler <- tibble(x = seq(2000, 26000, 2000),
                           y = 24.8,
                           yend = 25.3,
                           lab = paste0(seq(2, 26, 2), "kb")
    )
  }else if(reg == "E3"){
    genome_ruler <- tibble(x = seq(0, 26000, 2000),
                           y = 24.8,
                           yend = 25.3,
                           lab = paste0(seq(0, 26, 2), "kb")
    )
  }else{
    genome_ruler <- tibble(x = seq(0, 26266, 200),
                           y = 24.8,
                           yend = 25.3,
                           lab = paste0(seq(0, 26.266, 0.2), "kb")
    )
  }
  
  if(reg %in% c("E1")){
    offset_y <- 0.1
  }else if (reg %in% c("E2")){
    offset_y <- 0.15
  }else{offset_y <- 0.25}
 
  ljust <- ifelse(reg %in% c("E2", "E4"), 0, 1)
  rjust <- ifelse(reg %in% c("E2", "E4"), 1, 0)
  plot_full_trxptome(combined_gtf, reg, reg2) +
  # genome size marker
  geom_segment(data = genome_ruler, aes(x = x, xend = x, y = y, yend = yend),
               color = "#ffffff", linewidth = 0.5) +
  # genome size labels
  geom_richtext(data = genome_ruler, aes(x = x, y = (y - 0.1), label = lab),
                label.size = NA, label.padding = unit(0, "pt"), fill = NA) +
  # start codon positions
  geom_richtext(aes(x = sstcodon, y = ypos, label = glue("|")),
                size = 6, label.size = NA, label.padding = unit(0, "pt"),
                fill = NA, nudge_y = offset_y/2) +
  geom_richtext(aes(x = sstcodon, y = ypos, label = glue("<sup>SC ({sstcodon})</sup>")),
                size = 4.5, label.size = NA, label.padding = unit(0, "pt"),
                fill = NA, nudge_y = offset_y, hjust = ljust) +
  # stop codon positions
  geom_richtext(aes(x = stpcodon, y = ypos, label = glue("|")),
                size = 6, label.size = NA, label.padding = unit(0, "pt"),
                fill = NA, nudge_y = offset_y/2) +
    geom_richtext(aes(x = stpcodon, y = ypos, label = glue("<sup>STC ({stpcodon})</sup>")),
                  size = 4.5, label.size = NA, label.padding = unit(0, "pt"),
                  fill = NA, nudge_y = offset_y, hjust = rjust) +
    # secondary start codon positions
    geom_richtext(aes(x = secSSC, y = ypos, label = glue("|")),
                  size = 6, label.size = NA, label.padding = unit(0, "pt"),
                  fill = NA, nudge_y = -offset_y) +
    geom_richtext(aes(x = secSSC, y = ypos, label = glue("<sup>secSC ({secSSC})</sup>")),
                  size = 4.5, label.size = NA, label.padding = unit(0, "pt"),
                  fill = NA, nudge_y = -offset_y*1.5, hjust = ljust) +
    # secondary stop codon positions
    geom_richtext(aes(x = secSTC, y = ypos, label = glue("|")),
                  size = 6, label.size = NA, label.padding = unit(0, "pt"),
                  fill = NA, nudge_y = -offset_y) +
    geom_richtext(aes(x = secSTC, y = ypos, label = glue("<sup>secSTC ({secSTC})</sup>")),
                  size = 4.5, label.size = NA, label.padding = unit(0, "pt"),
                  fill = NA, nudge_y = -offset_y*1.5, hjust = rjust)

}

brkdwn_e1_trxtps <- brkdown_reg_plots("E1")

brkdwn_e2_trxtps <- brkdown_reg_plots(reg = "E2", reg2 = "IM")

brkdwn_e3_trxtps <- brkdown_reg_plots("E3")

brkdwn_e4_trxtps <- brkdown_reg_plots("E4")

# brkdwn_im_trxtps <- brkdown_reg_plots("IM")

brkdwn_mlp_trxtps <- brkdown_reg_plots("MLP")

# =========================================================================
# Table of transcripts with their primer and gel validation images
# for supplementary PRC methods section
# =========================================================================
trxpt_exons <- comp_spliced_gtf %>%
  mutate(full_trxpt = paste0(start_1, "-", end_1),
         exon1 = paste0(start_2, "-", end_2),
         exon2 = paste0(start_3, "-", end_3),
         exon3 = paste0(start_4, "-", end_4),
         exon4 = paste0(start_5, "-", end_5),
         exon5 = paste0(start_6, "-", end_6),
         exon6 = paste0(start_7, "-", end_7),
         exon7 = paste0(start_8, "-", end_8)) %>%
  select(trxpt_id, region, num_exons, full_trxpt, paste0("exon", c(1:7))) %>%
  mutate(exon3 = str_replace(exon3, "NA-NA", "-"),
         exon4 = str_replace(exon4, "NA-NA", "-"),
         exon5 = str_replace(exon5, "NA-NA", "-"),
         exon6 = str_replace(exon6, "NA-NA", "-"),
         exon7 = str_replace(exon7, "NA-NA", "-")
         # dplyr::across(c(base::paste0("exon", 3:7)), ~ str_replace(.x, "NA-NA", "-")),
         ) %>%
  arrange(region, trxpt_id)

wetlab_val <- tibble(trxpt_id = c(# E1
                                  paste0("TRXPT_", 1:4),
                                  # E2
                                  "TRXPT_21", "TRXPT_6", "TRXPT_7", "TRXPT_15", "TRXPT_31",
                                  # E3
                                  paste0("TRXPT_", c(22:27, 29)),
                                  # E4
                                  "TRXPT_28",
                                  # IM
                                  "TRXPT_5",
                                  paste0("TRXPT_", c(10:14, 16:20, 8:9))
                                  ),
                     forwardP = c(# E1
                                  "CCCggtaccTTCTGT\nTTGAATTGTGGGCGG",
                                  "CCCggtaccGAGGCCT\nGTTGGAATTGTTGC",
                                  "CCCggtacCATTTCCC\nGTACACGGTGTTG",
                                  "CCCggtaccGTCATCA\nCAACTGACCTTGTCGTC",
                                  # E2
                                  rep("CCCggtacCTGT\nTGCTGAGACTTCGGACC", 3), # E2 universal R,
                                  "CCCggtacCCTTTAAA\nATCAAGCCTATTGGTCTTGTAAC",
                                  "CCCggtacCTAGTGGC\nAGTGTTCGAAGATTCC",
                                  #E3
                                  rep("CCCggtacCTGA\nGGAGGTCGTAGACTCTGC", 4), #E3_univiersal F
                                  rep("CCCggtaccGTC\nCGAAGTCTCAGCAACAGATTC", 2),
                                  "CCCggtacCTGAGGAG\nGTCGTAGACTCTGC",
                                  # E4
                                  "CCCggtaccGGACAC\nGTGTTCGTTAGAGAACC",
                                  # IM
                                  "CCCggtaccTCTGGTGAGA\nTCTTCCAAACAGAAAG",
                                  # MLP
                                  rep("CCCggtaccGCTCATCATC\nCAGTTCTAAATTTCTCTCTGC", 5), #mlp_long_uni
                                  rep("CCCggtaccGGATCTC\nCAGATTCTGGTCTGTG", 3), #mlp_short_universal
                                  NA,
                                  "CCCggtaccGAGGATTTGA\nAGCCAATTATCCTTCAACG",
                                  rep("CCCggtaccGCTCATCAT\nCCAGTTCTAAATTTCTCTCTGC", 2)#mlp_long_uni
                                  ),
                     reverseP = c(# E1
                                  rep("CCCtctagaCGTCCA\nGTAGTCAGGAATTCTAGTG", 4),
                                  # E2A
                                  "CCCtctagaGAACCC\nAGATATTGGCTCCAAGG",
                                  # E2
                                  rep("CCCtctagaCATTGAATA\nGATAAGCGTAGCCAATCAGC", 2),
                                  "CCCtctagaGTGTCATT\nGTCTACGCTGTTGTAGTAG",
                                  "CCCtctagaCATTGCAGG\nTATGAATTGCGGAGTAG",
                                  # E3
                                  rep("CCCtctagaGCCA\nAGCTTGGTCAGGTGAC", 3), # E3 trxpt_B R
                                  "CCCtctagaGGTAGCACA\nTACTGTATTGCCTGAAGC",
                                  "CCCtctagaGCCAAG\nCTTGGTCAGGTGAC",
                                  "CCCtctagaTGCAAT\nGCTAATCCTCCTGCTG",
                                  "CCCtctagaGCCAAG\nCTTGGTCAGGTGAC",
                                  # E4
                                  "CCCtctagaCAGTG\nCAATCCGACGCTCTG",
                                  # IM
                                  "CCCtctagaCGCAA\nCCTGTAGGTCCGATTAC",
                                  #MLP
                                  "CCCtctagaCCTACTC\nTACGTCTCTTAGCAGC",
                                  "CCCtctagaGCTTCAG\nTATTAGCAGCTGCACAACC",
                                  "CCCtctagaTTTCC\nAGCTGAAGCCTGGAG",
                                  "CCCtctagaGCCAAG\nCTTGGTCAGGTGAC",
                                  "CCCtctagaGCTTCAGT\nATTAGCAGCTGCACAACC",
                                  "CCCtctagaGCCT\nGTCCAACAACCTGC",
                                  "CTCCCCATCTAGAC\nCTTTCATCTAACTG",
                                  "CCCtctagaGTTCTC\nCGTCTTCTACGTCGTG",
                                  NA,
                                  "CCCtctagaCTGCA\nGGCACAACAGGTG",
                                  "CCCtctagaCCTATC\nATCTGGCAATTCCGGTATG",
                                  "CCCtctagaCCTACT\nCTACGTCTCTTAGCAGC"
                                  ),
                     valid_status = c(rep("Validated", 3), "Not Validated",
                                      rep("Validated", 13+13)),
                     gel_image = c(# E1
                                   "wet_lab_validation/validation_gels/trxpt_1_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_2_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_3_gel.png",
                                   "Not Validated",
                                   # E2
                                   "wet_lab_validation/validation_gels/trxpt_21_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_6or7_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_6or7_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_15_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_31_gel.png",
                                   # E3
                                   "wet_lab_validation/validation_gels/trxpt_22or10j2_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_23or29_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_24or11j2_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_25_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_26_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_27_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_23or29_gel.png",
                                   # E4
                                   "wet_lab_validation/validation_gels/trxpt_28_gel.png",
                                   # IM
                                   "wet_lab_validation/validation_gels/trxpt_5_gel.png",
                                   #MLP
                                   "wet_lab_validation/validation_gels/trxpt_10or9_j1_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_14or11j1_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_12_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_13_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_14or11j1_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_16_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_17_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_18_gel.png",
                                   "N/A",
                                   "wet_lab_validation/validation_gels/trxpt_20_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_8_gel.png",
                                   "wet_lab_validation/validation_gels/trxpt_10or9_j1_gel.png"
                                   )
                     )

supp_pcr_meth_tab <- left_join(trxpt_exons, wetlab_val, by = "trxpt_id")

