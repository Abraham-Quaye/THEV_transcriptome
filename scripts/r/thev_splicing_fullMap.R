#!/usr/bin/env Rscript

library(magrittr)
library(ggtext)
library(glue)
library(rtracklayer)
library(ggbrace)
library(patchwork)
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
  rbind(tibble(transcript_id = "TCONS_00000029", strand = "+", start_1 = 18230, start_2 = 18230,
               start_3 = 20162, start_4 = NA, start_5 = NA, start_6 = NA, start_7 = NA, start_8 = NA,
               end_1 = 20732, end_2 = 18350, end_3 = 20732, end_4 = NA, end_5 = NA, end_6 = NA,
               end_7 = NA, end_8 = NA, region = "E3", trxpt_id = "TRXPT_29"))

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
  select(trxpt_id, t_name, num_exons) %>%
  rbind(tibble(trxpt_id = "TRXPT_29", t_name = "22K", num_exons = 2)) 

comp_spliced_gtf <- left_join(spliced_gtf, t_exp_levels,
                              by = "trxpt_id")


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
                            gene_name  %in% c("E3", "100K") ~ "E3",
                            gene_name == "IVa2" ~ "IM",
                            gene_name == "UXP" ~ "LONE",
                            TRUE ~ "MLP")) %>%
  arrange(strand, region)


combined_gtf <- plyr::rbind.fill(comp_spliced_gtf, predicted_orfs) %>% 
  arrange(strand, region) %>% as_tibble() %>%
  mutate(trxpt_id = ifelse(is.na(trxpt_id), gene_name, trxpt_id))



plot_full_trxptome <- function(combined_gtf, trxptome_part){
  
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
  ypos_pos_e3$ypos <- rep(c(36:45), length.out = nrow(ypos_pos_e3))
  
  
  # add y-axis positions to negative strand transcripts
  
  ypos_e4 <- combined_gtf %>% filter(region == "E4")
  ypos_e4$ypos <- rep(c(22, 23), length.out = nrow(ypos_e4))
  
  ypos_e2 <- combined_gtf %>% filter(region == "E2")
  ypos_e2$ypos <- rep(c(23:18), length.out = nrow(ypos_e2))
  
  ypos_im <- combined_gtf %>% filter(region == "IM")
  ypos_im$ypos <- rep(c(19, 20), length.out = nrow(ypos_im))
  
  ypos_uxp <- combined_gtf %>% filter(region == "LONE")
  ypos_uxp$ypos <- rep(c(23), length.out = nrow(ypos_uxp))
  
  
  trxpt_regs <- if(trxptome_part == "all"){
    c(unique(combined_gtf$region))  
  } else{c(trxptome_part)}
  
  # recreate full dataframe with y-axis positions
  combined_gtf <- rbind(ypos_pos_e1, ypos_pos_mlp, ypos_pos_e3,
                        ypos_e4, ypos_e2, ypos_im, ypos_uxp) %>%
    mutate(color = case_when(region == "E1" ~ "#ff0000",
                             region == "E2" ~ "#000000",
                             region == "E3" ~ "grey50",
                             region == "E4" ~ "#00ff00",
                             region == "IM" ~ "steelblue",
                             TRUE ~ "#0000ff")) %>%
    mutate(color = ifelse(!is.na(gene_name), "grey80", color),
           ypos = case_when(gene_name == "pVIII" ~ 29, 
                            gene_name == "Protease" ~ 43,
                            gene_name == "E3" ~ 42,
                            gene_name == "ORF7" ~ 26,
                            TRUE ~ ypos)) %>%
    as_tibble() %>%
    filter(region %in% trxpt_regs)
  
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
  
  # e3
  if(trxptome_part == "E3"){
    combined_gtf <- combined_gtf %>%
      mutate(ypos = rep(c(26:43), length.out = nrow(.)))
  }
  
  # ===================================================================
  # SINGLE EXON TRXPTS AND REGION LABELS FOR FULL TRANSCRIPTOME =======
  single_orfs <- combined_gtf %>% filter(is.na(start_2), is.na(end_2))
  
  genome_ruler <- tibble(x = seq(0, 26000, 2000),
                         y = 24.8,
                         yend = 25.3,
                         lab = paste0(seq(0,26, 2), "kb")
                         )
  
  trxpt_units <- tribble(~x, ~xend, ~y, ~labb, ~ylabb,
                         0, 2325, 17.5, "E1", 17,
                         24511, 26266, 17.5, "E4", 17,
                         2333, 3437, 17.5, "IM", 17,
                         3429, 8543, 17.5, "E2B", 17,
                         16972, 18186, 17.5, "E2A", 17
                         )
  e3_lab <- tribble(~x, ~xend, ~y, ~labb, ~ylabb,
                     18186, filter(combined_gtf, transcript_id == "TCONS_00000027") %>% pull(end_1), 44, "E3", 44.5
                    )
  
  mlp_lab <- tribble(~x, ~xend, ~y, ~labb, ~ylabb,
                     filter(combined_gtf, transcript_id == "TCONS_00000006") %>% pull(start_1), filter(combined_gtf, transcript_id == "TCONS_00000022") %>% pull(end_1), 44.5, "MLTU", 46.5)
  
  
  # ========================================================
  # plot transcripts 
  # ========================================================
  
  # Extract the column names that start with "start" or "end"
  start_cols <- names(combined_gtf)[str_detect(names(combined_gtf), "^start")]
  end_cols <- names(combined_gtf)[str_detect(names(combined_gtf), "^end")]
  
  # 
  x_axis_start <- min(combined_gtf[, start_cols[1]])
  if(trxptome_part != "all"){
    x_axis_start <- x_axis_start - (x_axis_start * 0.1)
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
    x_axis_end <- x_axis_end + (x_axis_end * 0.1)
    x_axis_end <- plyr::round_any(x_axis_end, 100, f = ceiling)
  } else{x_axis_end <- 26266}
  
  
  if(trxptome_part == "all"){
    y_lims <- c(NA, 48)
    } else{y_lims <- c(NA, NA)}
  
  splice_map <- combined_gtf %>%
    ggplot() +
    # line representing whole genome
    geom_segment(aes(x = x_axis_start, xend = x_axis_end, y = 25, yend = 25),
                 linewidth = 3.5, color = "black") +
   
    # plot full trxpts: start_pos to end_pos
    geom_segment(aes_string(x = start_cols[1], xend = end_cols[1],
                            y = "ypos", yend = "ypos"),
                 linetype = "dotted", color = combined_gtf$color) +
    
    # plot single orfs
    geom_segment(data = single_orfs, aes(x = start_1, xend = end_1, y = ypos, yend = ypos),
                 linewidth = 5, color = single_orfs$color) +
    
    # aesthetics
    scale_x_continuous(expand = c(0.01,0.01),
                       limits = c(x_axis_start, (x_axis_end + (x_axis_end * 0.005))),
                       breaks = seq(0, 26000, 2000), 
                       labels = paste0(seq(0,26, 2), "kb")) +
    scale_y_continuous(expand = c(0.01,0.01),
                       limits = y_lims) +
    theme(plot.margin = margin(rep(15, 4)),
          plot.background = element_rect(fill = "#ffffff"),
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
    splice_map <- splice_map +
      geom_segment(aes_string(x = start_cols[i], xend = end_cols[i], y = "ypos", yend = "ypos"),
                   color = combined_gtf$color, linewidth = 5)
  }
  
  # Add labels for full plot -----------------
  if(trxptome_part == "all"){
    splice_map <- splice_map +
      # genome size marker
      geom_segment(data = genome_ruler, aes(x = x, xend = x, y = y, yend = yend),
                   color = "#ffffff", linewidth = 0.5) +
      # genome size labels
      geom_richtext(data = genome_ruler, aes(x = x, y = (y - 0.5), label = lab),
                    label.size = NA, label.padding = unit(0, "pt"), fill = NA) +
      # transcription unit labels
      geom_segment(data = trxpt_units, aes(x = x, xend = xend, y = y, yend = y),
                   linewidth = 0.65, arrow = arrow(ends = "both", type = "closed",
                                                   angle = 25, length = unit(0.1, "inches"))) +
      geom_text(data = trxpt_units, aes(x = (x + xend)/2, y = ylabb, label = labb),
                fontface = "bold", size = 5) +
      
      # E3 transcription unit labels
      geom_text(data = e3_lab, aes(x = (x + xend)/2, y = ylabb, label = labb),
                fontface = "bold", size = 6) +
      
      geom_brace(aes(x = c(e3_lab$x, e3_lab$xend),
                     y = c(e3_lab$y, e3_lab$ylabb - 1)),
                 inherit.data = F) +
      
      # MLP transcription unit labels
      geom_text(data = mlp_lab, aes(x = (x + xend)/2, y = ylabb, label = labb),
                fontface = "bold", size = 6) +
      
      geom_brace(aes(x = c(mlp_lab$x, mlp_lab$xend),
                     y = c(mlp_lab$y, mlp_lab$ylabb - 1)),
                 inherit.data = F)
  }else{
    splice_map <- splice_map +
    # label each transcript
    geom_richtext(aes(x = (end_1), y = ypos, label = trxpt_id,
                      hjust = 0, fontface = "bold"), size = 3,
                  label.size = NA, label.padding = unit(0, "pt"),
                  fill = NA, label.margin = margin(l = 5),
                  show.legend = F)
  }
  
  if(trxptome_part == "all"){
    splice_map <- splice_map +
      plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(size = 22, face = "bold"))
  }
  return(splice_map)
}  

## Save final figure of full transcriptome
# ggsave(plot = plot_full_trxptome(combined_gtf, "all"),
#        filename = "results/r/figures/thev_spliced_map.png",
#        dpi = 500, width = 12, height = 8)



