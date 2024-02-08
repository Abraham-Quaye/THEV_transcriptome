#!/usr/bin/env Rscript


library(magrittr)
library(ggtext)
library(tidyverse)

# A function that loads in bed file, wrangles data to make plot-ready, and plots it

make_genomic_map <- function(bedfile){
  thev_genome <- read_tsv(bedfile,
                          col_names = FALSE,
                          col_types = "ciiciciicicc") %>% 
    set_colnames(c("chr", "start", "end", "gene_name", "score", "strand", "thickStart",
                   "thickEnd", "color", "blockCount", "blockSizes", "blockStarts")) %>% 
    mutate(exon1Size = as.numeric(str_replace(blockSizes, pattern = "(\\d+),\\d+", replacement = "\\1")),
           exon2Size = as.numeric(str_replace(blockSizes, pattern = "\\d+,(\\d+)", replacement = "\\1"))) %>%
    mutate(exon2Size = ifelse(exon2Size == exon1Size, NA_real_, exon2Size),
           exon1start = thickStart,
           exon1end = ifelse(blockCount > 1, (start + exon1Size), thickEnd),
           exon2start = thickEnd - exon2Size,
           exon2end = thickEnd,
           color = ifelse(strand == "+", "#FF0000", "skyblue"),
           gene_name = case_match(gene_name,
                                  "33K_spliced" ~ "33K",
                                  "pVIII gene" ~ "pVIII",
                                  .default = gene_name),
           ypos = c(rep(c(0.5, 1, 1.5),7), 0.5 , 1))
  
  spliced <- thev_genome %>% filter(blockCount > 1)
  
  ggplot(thev_genome) +
    # line representing whole genome
    geom_segment(aes(x = 0, xend = 26266, y = 0, yend = 0),
                 linewidth = 3.5, color = "black") +
    geom_segment(aes(x = thickStart, xend = thickEnd, y = ypos, yend = ypos),
                 color = thev_genome$color, linetype = "dashed",
                 ) +
    geom_segment(aes(x = exon1start, xend = exon1end, y = ypos, yend = ypos),
                 color = thev_genome$color, linewidth = 7) +
    geom_segment(data = spliced, aes(x = exon2start, xend = exon2end, y = ypos, yend = ypos),
                 color = spliced$color, linewidth = 7) +
    geom_text(aes(x = (thickStart + thickEnd)/2,
                   y = ypos,
                   label = gene_name),
               fontface = "bold", size = 3.5) +
    scale_x_continuous(expand = c(0, 0),
                       breaks = seq(0,26000,2000), 
                       labels = paste0(seq(0,26, 2), "kb")) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 1.6)) +
    theme(plot.margin = margin(rep(15, 4)),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(face = "bold", margin = margin(t = 5),
                                     size = 8.5),
          axis.ticks.length.x = unit(0.25, "cm"),
          axis.text.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.y = element_blank()
          )
}