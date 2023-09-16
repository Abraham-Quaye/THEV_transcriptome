#!/usr/bin/env Rscript

library(magrittr)
library(glue)
library(tidyverse)
library(rtracklayer)
library(patchwork)
library(ggtext)

# table of transcript ids for different genome regions here for transcription unit assignment

trxpt_ids_repo <- tribble(~region, ~`4hpi`, ~`12hpi`, ~`24hpi`, ~`72hpi`,
                     "e1", NA, c("12S1.1.2", "12S1.1.1"), c("24S1.1.1", "24S2.1.1",
                                                            "24S2.1.2", "24S3.1.2",
                                                            "24S3.1.1", "24S1.1.2"), c("72S1.1.1", "72S2.1.1", "72S2.1.2", "72S3.1.1"),
                     "e2", c("4S1.4.1",
                             "4S2.3.1"), NA, c("24S2.2.2", "24S3.2.2", "24S3.2.1",
                                               "24S1.2.1", "24S1.2.2", "24S2.2.1",
                                               "24S1.2.4", "24S3.2.5", "24S2.2.4"), c("72S1.12.1", "72S2.8.1", "72S3.10.1", "72S3.3.2"),
                     "e3", c("4S1.5.1"), c("12S3.6.2",
                                           "12S1.4.7",
                                           "12S1.4.8"), c("24S3.11.1", "24S1.11.1", "24S2.11.1",
                                                          "24S1.11.2", "24S3.11.2", "24S2.11.2"), c("72S1.13.4", "72S3.11.3", "72S1.13.1", "72S3.11.2", "72S1.13.2"),
                     "e4", NA, c("12S3.8.1"), NA, NA,
                     "im", NA, c("12S1.2.2"), NA, c("72S3.2.1")
                     ) %>% 
  pivot_longer(-region, names_to = "timepoint", values_to = "ids")
  
  duplicate_trxpts <- c("24S1.1.1", "24S1.1.2", "24S3.1.1", "72S3.1.1",
                        "72S1.1.1", "24S3.11.1", "24S1.11.1", "24S1.11.2",
                        "24S3.11.2", "72S1.13.1")
# "results/stringtie/transcripts_merged_72hrs.gtf"
# Import GTF and Wrangle for plotting ==========================================
#===============================================================================

graph_trxpts <- function(gtf){
  
  # detect timepoint from gtf name and pass onto next funtion 
  # to avail the appropriate trxpt ids
  catch_timepoint <- str_extract(gtf, "\\d+") %>% paste0(., "hpi")
  
  region_time_trxptIDs <- function(tp){
    e1_exon_names <- trxpt_ids_repo %>% filter(region == "e1", timepoint == tp) %>% pull(ids)
    e2_exon_names <- trxpt_ids_repo %>% filter(region == "e2", timepoint == tp) %>% pull(ids)
    e3_exon_names <- trxpt_ids_repo %>% filter(region == "e3", timepoint == tp) %>% pull(ids)
    e4_exon_names <- trxpt_ids_repo %>% filter(region == "e4", timepoint == tp) %>% pull(ids)
    im_exon_names <- trxpt_ids_repo %>% filter(region == "im", timepoint == tp) %>% pull(ids)

    return(tibble(e1 = e1_exon_names, e2 = e2_exon_names, e3 = e3_exon_names,
                  e4 = e4_exon_names, im = im_exon_names))
  }
  # this is a tibble containing trxpt_ids of all regions for the detected timepoint
  tp_trxpt_ids <- region_time_trxptIDs(catch_timepoint)
  
  # Now Wrangling Data for plotting ==========================================
  merged_gtf <- import(gtf) %>%
    as_tibble() %>%
    filter(!transcript_id %in% duplicate_trxpts) %>% 
    group_by(transcript_id, strand) %>%
    reframe(exon_count = n() - 1,
            start = list(start),
            end = list(end)) %>%
    unnest_wider(c(start, end), names_sep = "_") %>% 
    mutate(region = case_when(transcript_id %in% unlist(tp_trxpt_ids$e1) ~ "E1",
                              transcript_id %in% unlist(tp_trxpt_ids$e2) ~ "E2",
                              transcript_id %in% unlist(tp_trxpt_ids$e3) ~ "E3",
                              transcript_id %in% unlist(tp_trxpt_ids$e4) ~ "E4",
                              transcript_id %in% unlist(tp_trxpt_ids$im) ~ "IM",
                              TRUE ~ "MLP")) %>% 
    arrange(strand, region)
  
  # add y-axis positions to postive strand E1 transcripts
  ypos_vals <- c(26:45)
  ypos_pos_e1 <- merged_gtf %>% filter(strand == "+" & region == "E1")
  ypos_pos_e1$ypos <- rep(ypos_vals, length.out = nrow(ypos_pos_e1))
  
  # add y-axis positions to postive strand MLP transcripts
  ypos_pos_mlp <- merged_gtf %>% filter(strand == "+" & region == "MLP")
  ypos_pos_mlp$ypos <- rep(ypos_vals, length.out = nrow(ypos_pos_mlp))
  
  # add y-axis positions to postive strand E3 transcripts
  ypos_pos_e3 <- merged_gtf %>% filter(strand == "+" & region == "E3")
  ypos_pos_e3$ypos <- case_when(catch_timepoint == "4hpi" ~ rep(ypos_vals, length.out = nrow(ypos_pos_e3)),
                                catch_timepoint == "12hpi" ~ rep(c(27:32), length.out = nrow(ypos_pos_e3)),
                                catch_timepoint == "24hpi" ~ rep(ypos_vals, length.out = nrow(ypos_pos_e3)),
                                catch_timepoint == "72hpi" ~ rep(c(31:40), length.out = nrow(ypos_pos_e3))
                                )
  
  # add y-axis positions to negative strand transcripts
  if(catch_timepoint == "24hpi"){
    yneg_vals <- c(22:24)
    }else{yneg_vals <- c(22)}
  
  ypos_neg_gtf <- merged_gtf %>% filter(strand == "-")
  ypos_neg_gtf$ypos <- rep(yneg_vals, length.out = nrow(ypos_neg_gtf))
  
  # recreate full dataframe with y-axis positions
  merged_gtf <- rbind(ypos_pos_e1, ypos_pos_mlp, ypos_pos_e3, ypos_neg_gtf) %>% 
    mutate(color = case_when(region == "E1" ~ "#ff0000",
                             region == "E2" ~ "#000000",
                             region == "E3" ~ "grey50",
                             region == "E4" ~ "#00ff00",
                             region == "IM" ~ "steelblue",
                             TRUE ~ "#0000ff"))
  
  # data for plotting line representing full thev genome
  genome_ruler <- tibble(x = seq(0,26000,2000),
                         y = 24.95,
                         yend = 25.05,
                         lab = paste0(seq(0,26, 2), "kb")
  )
  
  #### -----------------------------------------------------------------
  # Extract the column names that start with "start" or "end"
  # Plot All exons until complete
  # ====================================================================
    ylimits <- case_when(catch_timepoint == "4hpi" ~ c(21, 26.5),
                         catch_timepoint == "12hpi" ~ c(21, 30),
                         catch_timepoint == "24hpi" ~ c(21, 29),
                         catch_timepoint == "72hpi" ~ c(21, 35)
    )
  
  
    plot_tibble <- function(df) {
      # Extract the column names that start with "start" or "end"
      start_cols <- names(df)[str_detect(names(df), "^start")]
      end_cols <- names(df)[str_detect(names(df), "^end")]
      
      # Create the initial plot with the first pair of start and end columns
      plott <- ggplot(df) +
        # line representing whole genome
        geom_segment(aes(x = 0, xend = 26266, y = 25, yend = 25),
                     linewidth = 1, color = "black") +
        # genome size marker
        geom_segment(data = genome_ruler, aes(x = x, xend = x, y = y, yend = yend),
                     color = "#ffffff", linewidth = 0.5) +
        
        # genome size labels
        geom_richtext(data = genome_ruler, aes(x = x, y = (y - 0.3), label = lab),
                      label.size = NA, label.padding = unit(0, "pt"), fill = NA) +
        
        # plot full trxpt outline as dotted line: start_pos to end_pos
        geom_segment(aes_string(x = start_cols[1], xend = end_cols[1], y = "ypos", yend = "ypos"),
                     linetype = "dotted", color = df$color) +
        # identify gtf being plotted
        labs(title = glue("Transcript map at {catch_timepoint}")) +
        # aesthetics
        scale_x_continuous(expand = c(0.01,0.01),
                           limits = c(0, 26400),
                           breaks = seq(0,26000,2000),
                           labels = paste0(seq(0,26, 2), "kb")) +
        ## to plot the plot independently, use limits = c(15, 50) for best results
        scale_y_continuous(expand = c(0.01,0.01),
                           limits = ylimits) +
        theme(plot.margin = margin(rep(5, 4)),
              plot.background = element_rect(fill = "#ffffff", color = "grey"),
              panel.background = element_rect(fill = "#ffffff"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.title = element_blank(),
              axis.text= element_blank(),
              axis.line.x = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(hjust = 0.1, colour = "#000000",
                                        size = 10, face = "bold"))
      
      # Loop through the remaining start and end columns and 
      # add a geom_segment layer for each pair
      for (i in seq_along(start_cols)[-1]) {
        plott <- plott +
          geom_segment(aes_string(x = start_cols[i], xend = end_cols[i], y = "ypos", yend = "ypos"),
                       color = df$color, linewidth = 5)
      }
      # Return the final plot
      return(plott)
    }
 f_plot <- plot_tibble(merged_gtf)
 return(f_plot)
}

# load gtf files from source for plotting ======================================

tp_trxptome <- list.files("results/stringtie",
                          pattern = "transcripts_merged_\\d{1,2}hrs\\.gtf",
                          full.names = T)

# Call function to plot full transcriptome at each time point and save using a for loop
# ==============================================================================

# for(gtf_file in tp_trxptome){
#   catch_timepoint <- str_extract(gtf_file, "\\d+") %>% paste0(., "hpi")
#   graph_trxpts(gtf_file)
#   ggsave(filename = glue("results/r/figures/thev_spliced_map_{catch_timepoint}.png"),
#          dpi = 500, width = 12, height = 8)
# }

p12 <- graph_trxpts(tp_trxptome[1])
p24 <- graph_trxpts(tp_trxptome[2])
p4 <- graph_trxpts(tp_trxptome[3])
p72 <- graph_trxpts(tp_trxptome[4])

all_plots <- (p4/ p12 / p24 / p72) +
  plot_layout(tag_level = "new",
              heights = c(2, 3, 3, 4)) +
  plot_annotation(tag_levels = "1", tag_prefix = "B") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(plot = all_plots,
       filename = "results/r/figures/thev_patched_timepoints_spliced_map.png",
       dpi = 500, width = 12, height = 10.5)


