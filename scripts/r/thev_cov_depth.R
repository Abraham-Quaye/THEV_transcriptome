#!/usr/bin/env Rscript

library(magrittr)
library(glue)
library(patchwork)
library(ggsci)
library(ggtext)
library(grDevices)
library(tidyverse)

source("scripts/r/thev_genomic_map.R")
# load in depth files
genome <- make_genomic_map("raw_files/annotations/THEVannotated_genesOnly.bed")
### find all the depth files
all_depth_files <- list.files("results/hisat2/coverage",
                              pattern = "thev_\\d+[a-z]+\\.txt",
                              full.names = TRUE) %>% 
  setNames(c("12hpi", "24hpi", "4hpi", "72hpi"))

# read depth_files as one master tibble
all_depths <- map_dfr(all_depth_files, read_tsv,
                      col_names = FALSE,
                      show_col_types = FALSE,
                      comment = "#",
                      .id = "timepoint") %>% 
  mutate(timepoint = factor(timepoint,
                            levels = c("72hpi", "24hpi", "12hpi", "4hpi")),
         color = c(rep("grey50",26266),
                   rep("grey40",26266),
                   rep("grey",26266),
                   rep("grey30",26266)
                   ),
         titles = c(rep("12hrs Post-infection", 26266),
                    rep("24hrs Post-infection", 26266),
                    rep("4hrs Post-infection", 26266),
                    rep("72hrs Post-infection", 26266)
                    )) %>% 
  set_colnames(c("timepoint", "genome", "position",
                 "depth", "color", "titles"))

## --------------------
# function to prepare each timepoint for plotting on a log scale
log_plot_times <- function(tp){
  tp_data <- all_depths %>%
  select(timepoint, position, depth) %>%
  filter(timepoint == tp) %>%
  mutate(depth = depth + 1) %>%
  mutate(log_depth = log10(depth))
  
  return(tp_data)
}
log_data_tps <- map(c("72hpi", "24hpi", "12hpi", "4hpi"), log_plot_times)


comp_all <- ggplot() +
  geom_area(data = log_data_tps[[1]],
            aes(position, log_depth, group = 1, fill = "72hpi"),
            key_glyph = "crossbar") +
  geom_area(data = log_data_tps[[2]],
            aes(position, log_depth, group = 2, fill = "24hpi"),
            key_glyph = "crossbar") +
  geom_area(data = log_data_tps[[3]],
            aes(position, log_depth, group = 3, fill = "12hpi"),
            key_glyph = "crossbar") +
  geom_area(data = log_data_tps[[4]],
            aes(position, log_depth, group = 4, fill = "4hpi"),
            key_glyph = "crossbar") +
  scale_fill_manual(values = c("#5050ff", "#ce3d32", "#749b58", "#f0e685"),
                    breaks = c("72hpi", "24hpi", "12hpi", "4hpi"),
                    guide = guide_legend(keywidth = unit(1, "cm"))) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(0, 26000, 1000),
                     labels = glue("{seq(0, 26, 1)}kb")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Mapping Depth of RNA-seq Reads Over THEV Genome",
       x = element_blank(),
       y = "Mapping Depth (Log10)",
       fill = element_blank(),
       color = element_blank()) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 16,
                                   face = "bold",
                                   color = "#000000",
                                   margin = margin(l = 10, r = 5)),
        axis.text.x = element_text(size = 12,
                                   face = "bold",
                                   color = "#000000",
                                   margin = margin(b = 10, t = 5)),
        plot.title = element_text(size = 28, 
                                  face = "bold", 
                                  hjust = 0.5,
                                  margin = margin(t = 10)),
        axis.title.y = element_text(size = 22, face = "bold",
                                    margin = margin(l = 10)),
        legend.justification = c(0, 0),
        legend.position = c(0, 0.92),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.text = element_text(size = 14,
                                   margin = margin(r = 10, l = 0, t = 0, b = 0),
                                   face = "bold"),
        legend.key.height = unit(1, "cm"),
        legend.text.align = 0)
