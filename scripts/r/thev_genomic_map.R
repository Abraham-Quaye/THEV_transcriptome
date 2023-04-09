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
    mutate(exon1Size = str_replace(blockSizes, pattern = "(\\d+),\\d+", replacement = "\\1"),
           exon1Size = as.numeric(exon1Size),
           exon2Size = str_replace(blockSizes, pattern = "\\d+,(\\d+)", replacement = "\\1"),
           exon2Size  = as.numeric(exon2Size),
           exon2Size = ifelse(exon2Size == exon1Size, NA_real_, exon2Size),
           exon1start = thickStart,
           exon1end = ifelse(blockCount > 1, (start + exon1Size), thickEnd),
           exon2start = thickEnd - exon2Size,
           exon2end = thickEnd,
           color = ifelse(strand == "+", "#B897A4", "#8B62D1"),
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



## position the plus strand genes above genome line and minus strand genes below
make_genomic_map2 <- function(bedfile){
  thev_genome2 <- read_tsv(bedfile,
                          col_names = FALSE,
                          col_types = "ciiciciicicc") %>% 
    set_colnames(c("chr", "start", "end", "gene_name", "score", "strand", "thickStart",
                   "thickEnd", "color", "blockCount", "blockSizes", "blockStarts")) %>% 
    mutate(exon1Size = str_replace(blockSizes, pattern = "(\\d+),\\d+", replacement = "\\1"),
           exon1Size = as.numeric(exon1Size),
           exon2Size = str_replace(blockSizes, pattern = "\\d+,(\\d+)", replacement = "\\1"),
           exon2Size  = as.numeric(exon2Size),
           exon2Size = ifelse(exon2Size == exon1Size, NA_real_, exon2Size),
           exon1start = thickStart,
           exon1end = ifelse(blockCount > 1, (start + exon1Size), thickEnd),
           exon2start = thickEnd - exon2Size,
           exon2end = thickEnd,
           color = ifelse(strand == "+", "#B897A4", "#8B62D1"),
           gene_name = case_match(gene_name,
                                  "33K_spliced" ~ "33K",
                                  "pVIII gene" ~ "pVIII",
                                  .default = gene_name),
           ypos = ifelse(strand == "-", c(2, 3), c(6, 7)),
           ypos = ifelse(gene_name == "22K", 8, ypos))
  
  spliced2 <- thev_genome2 %>% filter(blockCount > 1)
  
  genome_ruler <- tibble(x = seq(0,26000,2000),
                         xend = seq(0,26000,2000),
                         y = 4.7,
                         yend = 5.3,
                         lab = paste0(seq(0,26, 2), "kb"))
  
  
  ggplot(thev_genome2) +
    # line representing whole genome
    geom_segment(aes(x = 0, xend = 26266, y = 5, yend = 5),
                 linewidth = 3.5, color = "black") +
    # genome size marker
    geom_segment(data = genome_ruler, aes(x = x, xend = xend, y = y, yend = yend),
                 color = "#ffffff", linewidth = 0.5) +
    # genome size labels
    geom_richtext(data = genome_ruler, aes(x = x, y = (y - 0.5), label = lab),
                  label.size = NA, label.padding = unit(0, "pt")) +
    # plot genes
    geom_segment(aes(x = thickStart, xend = thickEnd, y = ypos, yend = ypos),
                 color = thev_genome2$color, linetype = "dashed") +
    geom_segment(aes(x = exon1start, xend = exon1end, y = ypos, yend = ypos),
                 color = thev_genome2$color, linewidth = 7) +
    geom_segment(data = spliced2, aes(x = exon2start, xend = exon2end, y = ypos, yend = ypos),
                 color = spliced2$color, linewidth = 7) +
    geom_text(aes(x = (thickStart + thickEnd)/2, y = ypos, label = gene_name),
              fontface = "bold", size = 3.5) +
    scale_x_continuous(expand = c(0.01,0.01),
                       limits = c(0, 26400),
                       breaks = seq(0,26000,2000), 
                       labels = paste0(seq(0,26, 2), "kb")) +
    scale_y_continuous(expand = c(0.01,0.01),
                       limits = c(0, 50)) +
    theme(plot.margin = margin(rep(15, 4)),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.title = element_blank(),
          axis.text= element_blank(),
          axis.line.x = element_blank(),
          axis.ticks = element_blank(),
    )
}
