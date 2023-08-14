
library(magrittr)
library(ggtext)
library(glue)
library(rtracklayer)
library(ggbrace)
library(patchwork)
library(plyr)
library(tidyverse)

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
                            gene_name %in% c("ORF8", "ORF7") ~ "E4",
                            gene_name %in% c("DBP", "pTP", "AdPOL") ~ "E2",
                            gene_name  %in% c("E3", "100K") ~ "E3",
                            gene_name == "IVa2" ~ "IM",
                            gene_name == "UXP" ~ "LONE",
                            TRUE ~ "MLP")) %>%
  arrange(strand, region)


combined_gtf <- plyr::rbind.fill(spliced_gtf, predicted_orfs) %>% 
  arrange(strand, region) %>% as_tibble()






plot_full_trxptome <- function(combined_gtf){
  
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
  ypos_e4$ypos <- rep(c(23, 24), length.out = nrow(ypos_e4))
  
  ypos_e2 <- combined_gtf %>% filter(region == "E2")
  ypos_e2$ypos <- rep(c(23:18), length.out = nrow(ypos_e2))
  
  ypos_im <- combined_gtf %>% filter(region == "IM")
  ypos_im$ypos <- rep(c(19, 20), length.out = nrow(ypos_im))
  
  ypos_uxp <- combined_gtf %>% filter(region == "LONE")
  ypos_uxp$ypos <- rep(c(23), length.out = nrow(ypos_uxp))
  
  
  # recreate full dataframe with y-axis positions
  combined_gtf <- rbind(ypos_pos_e1, ypos_pos_mlp, ypos_pos_e3,
                        ypos_e4, ypos_e2, ypos_im, ypos_uxp) %>% 
    mutate(color = case_when(region == "E1" ~ "#ff0000",
                             region == "E2" ~ "#000000",
                             region == "E3" ~ "grey50",
                             region == "E4" ~ "#00ff00",
                             region == "IM" ~ "steelblue",
                             TRUE ~ "#0000ff")) %>%
    mutate(color = ifelse(!is.na(gene_name), "grey90", color),
           ypos = case_when(gene_name == "pVIII" ~ 29, 
                            gene_name == "Protease" ~ 43,
                            gene_name == "E3" ~ 42,
                            TRUE ~ ypos)) %>%
    as_tibble()
  
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
                     18186, filter(combined_gtf, transcript_id == "TCONS_00000027") %>% pull(end_1), 43, "E3", 43.5
                    )
  
  mlp_lab <- tribble(~x, ~xend, ~y, ~labb, ~ylabb,
                     filter(combined_gtf, transcript_id == "TCONS_00000006") %>% pull(start_1), filter(combined_gtf, transcript_id == "TCONS_00000022") %>% pull(end_1), 44.5, "MLP", 46.5)
  
  
  # plot transcripts =================================================
  # ========================================================
  
  splice_map <- combined_gtf %>% 
    ggplot() +
    # line representing whole genome
    geom_segment(aes(x = 0, xend = 26266, y = 25, yend = 25),
                 linewidth = 3.5, color = "black") +
    # genome size marker
    geom_segment(data = genome_ruler, aes(x = x, xend = x, y = y, yend = yend),
                 color = "#ffffff", linewidth = 0.5) +
    # genome size labels
    geom_richtext(data = genome_ruler, aes(x = x, y = (y - 0.5), label = lab),
                  label.size = NA, label.padding = unit(0, "pt"), fill = NA) +
    
    # plot full trxpts: start_pos to end_pos
    geom_segment(aes(x = start_1, xend = end_1, y = ypos, yend = ypos),
                 linetype = "dotted", color = combined_gtf$color) +
    
    # plot exons
    geom_segment(aes(x = start_2, xend = end_2, y = ypos, yend = ypos),
                 linewidth = 5, color = combined_gtf$color) +
    
    geom_segment(aes(x = start_3, xend = end_3, y = ypos, yend = ypos),
                 linewidth = 5, color = combined_gtf$color) +
    
    geom_segment(aes(x = start_4, xend = end_4, y = ypos, yend = ypos),
                 linewidth = 5, color = combined_gtf$color) +
    
    geom_segment(aes(x = start_5, xend = end_5, y = ypos, yend = ypos),
                 linewidth = 5, color = combined_gtf$color) +
    
    geom_segment(aes(x = start_6, xend = end_6, y = ypos, yend = ypos),
                 linewidth = 5, color = combined_gtf$color) +
    
    geom_segment(aes(x = start_7, xend = end_7, y = ypos, yend = ypos),
                 linewidth = 5, color = combined_gtf$color) +
    
    geom_segment(aes(x = start_8, xend = end_8, y = ypos, yend = ypos),
                 linewidth = 5, color = combined_gtf$color) +
    
    # plot single orfs
    geom_segment(data = single_orfs, aes(x = start_1, xend = end_1, y = ypos, yend = ypos),
                 linewidth = 5, color = single_orfs$color) +
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
               inherit.data = F) +
    
    # aesthetics
    scale_x_continuous(expand = c(0.01,0.01),
                       limits = c(0, 26400),
                       breaks = seq(0,26000,2000), 
                       labels = paste0(seq(0,26, 2), "kb")) +
    scale_y_continuous(expand = c(0.01,0.01),
                       limits = c(NA, 48)) +
    theme(plot.margin = margin(rep(15, 4)),
          plot.background = element_rect(fill = "#ffffff"),
          panel.background = element_rect(fill = c("#FFF0F0")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.title = element_blank(),
          axis.text= element_blank(),
          axis.line.x = element_blank(),
          axis.ticks = element_blank(),
    )
  
  splice_map <- splice_map +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 22, face = "bold"))
  
  return(splice_map)
}  


ggsave(plot = plot_full_trxptome(combined_gtf),
       filename = "results/r/figures/thev_spliced_map.png",
       dpi = 1000, width = 12, height = 8)


