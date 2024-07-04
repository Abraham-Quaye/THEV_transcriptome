#!/usr/bin/env Rscript


library(magrittr)
library(ggtext)
library(tidyverse)

## position the plus strand genes above genome line and minus strand genes below
make_genomic_map2 <- function(bedfile){
  genus_specific <- c("ORF1", "Hyd", "E3", "ORF7", "ORF8")
  
  thev_genome2 <- read_tsv(bedfile,
                          col_names = FALSE,
                          col_types = "ciiciciicicc") %>% 
    set_colnames(c("chr", "start", "end", "gene_name", "score", "strand", "thickStart",
                   "thickEnd", "color", "blockCount", "blockSizes", "blockStarts")) %>% 
    mutate(exon1Size = as.numeric(str_replace(blockSizes, pattern = "(\\d+),\\d+", replacement = "\\1")),
           exon2Size = as.numeric(str_replace(blockSizes, pattern = "\\d+,(\\d+)", replacement = "\\1")),
           exon2Size = ifelse(exon2Size == exon1Size, NA_real_, exon2Size),
           exon1start = thickStart,
           exon1end = ifelse(blockCount > 1, (start + exon1Size), thickEnd),
           exon2start = thickEnd - exon2Size,
           exon2end = thickEnd,
           color = case_when(gene_name %in% genus_specific ~ "red",
                             gene_name == "U exon" ~ "skyblue",
                             gene_name == "22K" ~ "green",
                             TRUE ~ "#6D58F5"),
           gene_name = case_match(gene_name,
                                  "33K_spliced" ~ "33K",
                                  "pVIII gene" ~ "pVIII",
                                  .default = gene_name),
           ypos = ifelse(strand == "-", c(24, 24.5), c(25.5, 26))) %>%
    mutate(ypos = case_when(gene_name %in% c("22K", "pVI") ~ 26.5, 
                            gene_name == "pVIII" ~ 27,
                            TRUE ~ ypos),
           offpos = ifelse(strand == "-", -0.2, 0.2))
  
  spliced2 <- thev_genome2 %>% filter(blockCount > 1)
  
  genome_ruler <- tibble(x = seq(0,26000,2000),
                         y = 24.9,
                         yend = 25.1,
                         lab = paste0(seq(0,26, 2), "kb"))
  
  shade_regions <- tribble(~xmin, ~xmax, ~ymin, ~ymax,
                           0, thev_genome2$end[2], 23, 27,
                           thev_genome2$start[19], thev_genome2$end[19], 23, 27,
                           thev_genome2$start[22], 26266, 23, 27
                           )
  
  trxpt_units <- tribble(~x, ~xend, ~y, ~labb, ~ylabb,
                         0, thev_genome2$end[2], 23, "E1", 22.85,
                         thev_genome2$start[19], thev_genome2$end[19], 23, "E3", 22.85,
                         thev_genome2$start[22], 26266, 23, "E4", 22.85,
                         thev_genome2$start[3], thev_genome2$end[3], 23, "IM", 22.85,
                         thev_genome2$start[4], thev_genome2$exon1end[5], 23, "E2B", 22.85,
                         thev_genome2$start[14], thev_genome2$end[14], 23, "E2A", 22.85
                         )

  ggplot(thev_genome2) +
    geom_rect(data = shade_regions,
              aes(NULL, NULL, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              alpha = 0.15) +
    # line representing whole genome
    annotate(geom = "segment", x = 0, xend = 26266, y = 25, yend = 25,
                 linewidth = 3.5, color = "black") +
    # genome size marker
    geom_segment(data = genome_ruler, aes(x = x, xend = x, y = y, yend = yend),
                 color = "#ffffff", linewidth = 0.5) +
    # genome size labels
    geom_richtext(data = genome_ruler, aes(x = x, y = (y - 0.2), label = lab),
                  label.size = NA, label.padding = unit(0, "pt"), fill = NA) +
    # plot genes
    geom_segment(aes(x = thickStart, xend = thickEnd, y = ypos, yend = ypos),
                 color = thev_genome2$color, linetype = "dashed") +
    geom_segment(aes(x = exon1start, xend = exon1end, y = ypos, yend = ypos),
                 color = thev_genome2$color, linewidth = 7) +
    geom_segment(data = spliced2, aes(x = exon2start, xend = exon2end, y = ypos, yend = ypos),
                 color = spliced2$color, linewidth = 7) +
    # gene names
    geom_text(aes(x = (thickStart + thickEnd)/2, y = ypos, label = gene_name),
              fontface = "bold", size = 5, nudge_y = thev_genome2$offpos) +
    # transcription unit labels
    geom_segment(data = trxpt_units, aes(x = x, xend = xend, y = y, yend = y),
                 linewidth = 1, arrow = arrow(ends = "both", type = "closed",
                                              angle = 25, length = unit(0.1, "inches"))) +
    geom_text(data = trxpt_units, aes(x = (x + xend)/2, y = ylabb, label = labb),
              fontface = "bold", size = 5) +
    scale_x_continuous(expand = c(0.01,0.01),
                       limits = c(0, 26400),
                       breaks = seq(0,26000,2000), 
                       labels = paste0(seq(0,26, 2), "kb")) +
    scale_y_continuous(expand = c(0.01,0.01),
                       limits = c(22, 27.5)) +
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
}

predicted_genemap <- make_genomic_map2("raw_files/annotations/THEVannotated_genesOnly.bed")

ggsave(plot = predicted_genemap, filename = "results/r/figures/figure1.png",
       dpi = 300, width = 12, height = 7)
