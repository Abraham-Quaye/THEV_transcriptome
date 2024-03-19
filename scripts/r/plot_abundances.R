#!/usr/bin/env Rscript

library(magrittr)
library(scales)
library(tidyverse)
library(plyr)
library(rtracklayer)
library(glue)
library(ggsci)
library(patchwork)

## load in analyzed data
source("scripts/r/abundance_analyses.R")


## visualize ---------------------------
plot_time_genexpression <- function(count_ss){
  up_lim = max(count_ss$reg_junc_percent) %>% plyr::round_any(., 10, f = ceiling)
  count_ss %>%
    filter(timepoint != "4hpi") %>%
    ggplot(aes(timepoint, reg_junc_percent, group = region, color = region)) +
    geom_line(linewidth = 3) +
    geom_point(size = 5) +
    scale_colour_jco() +
    scale_y_continuous(expand = c(0.01, 0.01),
                       breaks = c(seq(0, up_lim, 10)),
                       labels = scales::label_percent(scale = 1)) +
    scale_x_discrete(expand = c(0.005, 0.01)) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(plot.margin = margin(rep(30, 4)),
          panel.grid.major.y = element_line(linewidth = 0.4, color = "grey",
                                            linetype = "dashed"),
          axis.title.y = element_text(size = 28, face = "bold", margin = margin(r = 15)),
          axis.text.x = element_text(size = 26, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"),
          legend.justification = c(0, 0),
          legend.position = c(0.05, 0.9),
          legend.key.width = unit(2, "cm"),
          legend.title = element_blank(),
          legend.text = element_text(size = 20, face = "bold", colour = "black"),
          legend.text.align = 0) +
    guides(color = guide_legend(label.position = "top",
                                label.hjust = 0.5,
                                label.vjust = 1,
                                direction = "horizontal",
                                nrow = 1,
                                ))
}


# plot and save abundances of juncs found in final trxptome
trxptome_juncs <- plot_time_genexpression(bulk_count_trxptome_ss) +
  labs(x = element_blank(),
       y = "Relative Abundances of Junctions in Transcriptome")


all_juncs <- plot_time_genexpression(plot_all_junc_abunds) +
  labs(x = element_blank(),
       y = "Relative Abundances of All Junctions")

# -------------------
## splice donor and acceptor frequencies

# acceptors and donors in trxptome
bulk_trxptome_ss_seq <- count_trxptome_ss(unq_bulk_juncs) %>%
  distinct(junc_ss, junc_ts, .keep_all = T) %>%
  dplyr::group_by(splice_site) %>%
  reframe(freq = n()) %>% 
  mutate(sum_junc_count = sum(freq),
         percent_abund = (freq / sum_junc_count) * 100)


# All acceptors and donors
all_ss_seq_ready <- bulk_junc_stats %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("4hpi", "12hpi",
                                       "24hpi", "72hpi"))) %>%
  group_by(timepoint, splice_site) %>%
  reframe(tally = n(),
          total_reads = sum(read_count)) %>%
  split(.$timepoint) %>%
  map(mutate, tot_tally_tp = sum(tally)) %>%
  do.call("rbind", .) %>%
  mutate(percent_abund = (tally / tot_tally_tp) * 100)

# function to make splice site plots
plot_tp_ss <- function(ss_data_tp = all_ss_seq_ready, tp){
  # filter out needed tp
  tp_only <- ss_data_tp %>% 
    filter(timepoint == tp)
  
  # get timepoint
  timpnt <- tp_only %>% pull(timepoint) %>% unique(.)
  
  # set plot parameters
  if(timpnt == "72hpi"){
    pnt_size <- 3
    fnt_size <- 3
    x_ax_size <- 6
    p_col <- "red"
    ylimm <- 40
    xpand <- 0.01
    push_y <- 0.1 * ylimm
  }else if(timpnt %in% c("24hpi")){
    pnt_size <- 5
    fnt_size <- 3.5
    x_ax_size <- 9
    p_col <- "blue"
    ylimm <- 80
    xpand <- 0.025
    push_y <- 0.1 * ylimm
  }else{
    pnt_size <- 7
    fnt_size <- 3
    x_ax_size <- 12
    p_col <- base::sample(c("grey30", "orange"), 1)
    ylimm <- 105
    xpand <- 0.05
    push_y <- 0.075 * ylimm
  }
  
  # make plot
  tp_only %>%
  ggplot(aes(splice_site, percent_abund,
                  group = timepoint, color = timepoint)) +
  geom_point(show.legend = F, size = pnt_size) +
  geom_segment(aes(x = splice_site, xend = splice_site,
                   y = 0, yend = percent_abund),
               linewidth = 2,
               show.legend = F) +
  geom_text(aes(label = glue("{round(percent_abund, 1)}%")),
            nudge_y = push_y, size = fnt_size, angle = 90,
            fontface = "bold", color = "#000000") +
  labs(title = base::paste("Unique Splice Sites Detected at", tp),
       x = "Splice Site Donor-Acceptor",
       y = "Percentage Frequency") +
  scale_color_manual(values = p_col) +
  scale_y_continuous(expand = c(0.01,0.01),
                     limits = c(0, ylimm),
                     labels = scales::label_percent(scale = 1)) +
  scale_x_discrete(expand = c(xpand, xpand)) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(plot.margin = margin(rep(30, 4)),
        panel.grid.major.y = element_line(linewidth = 0.4, color = "grey",
                                          linetype = "dashed"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 22, face = "bold",
                                  margin = margin(r = 15, t = 15)),
        axis.text.x = element_text(size = x_ax_size, face = "bold",
                                   angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(color = "grey40"),
        axis.ticks.length = unit(0, "pt"),
        panel.spacing = unit(1, "cm")
  )
}

timpoints <- paste0(c(4, 12, 24, 72), "hpi")
ss_all_ploted <- list()
for(p in seq_along(timpoints)){
  tpp <- timpoints[p]
  ss_all_ploted[[tpp]] <- plot_tp_ss(tp = tpp)
}

temp_ss_all <- (ss_all_ploted$`4hpi` | ss_all_ploted$`12hpi` | ss_all_ploted$`24hpi`) +
  plot_layout(widths = c(1,1.2,2))

save_ss_all <- temp_ss_all/ss_all_ploted$`72hpi` +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(plot = save_ss_all, filename = "results/r/figures/figure_5a_d.png",
       width = 25, height = 15, dpi = 350)

# ====================================================================
# PLOT TRANSCRIPT ABUNDANCES FROM BALLGOWN
# ====================================================================

# plot transcript expression levels per individual
trxpt_fpkms <- t_expr_each %>%
  filter(timepoint != "4h.p.i") %>%
  ggplot(aes(trxpt_id, fpkm_trxpt_percent, group = timepoint, fill = timepoint)) +
  geom_col(position = position_dodge(0.5), width = 0.5) +
  scale_fill_jco() +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = scales::label_percent(scale = 1)) +
  scale_x_discrete(expand = c(0.005, 0.01)) +
  coord_cartesian(clip = "off") +
  labs(x = element_blank(),
       y = "Proportion FPKM/Transcript (%)") +
  theme_classic() +
  theme(plot.margin = margin(t = 20, l = 20, r = 20, b = 5),
        panel.grid.major.y = element_line(linewidth = 0.2, color = "grey",
                                          linetype = "dashed"),
        axis.title.y = element_text(size = 28, face = "bold", margin = margin(r = 15)),
        axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 28, face = "bold", colour = "#000000"),
        legend.justification = c(0, 0),
        legend.position = c(0.02, 0.6),
        legend.key.size = unit(1.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, face = "bold", colour = "black"),
        legend.text.align = 0)

# plot transcript expression levels per region
reg_fpkms <- t_exp_lev_byregion %>%
  filter(timepoint != "4h.p.i") %>%
  ggplot(aes(timepoint, fpkm_reg_percent, group = region, color = region)) +
  geom_line(linewidth = 3) +
  geom_point(size = 5) +
  scale_colour_jco() +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = scales::label_percent(scale = 1)) +
  scale_x_discrete(expand = c(0.005, 0.01)) +
  coord_cartesian(clip = "off") +
  labs(x = element_blank(),
       y = "Proportion FPKM/Region (%)") +
  theme_classic() +
  theme(plot.margin = margin(rep(30, 4)),
        panel.grid.major.y = element_line(linewidth = 0.4, color = "grey",
                                          linetype = "dashed"),
        axis.title.y = element_text(size = 28, face = "bold", margin = margin(r = 15)),
        axis.text.x = element_text(size = 28, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        legend.justification = c(0, 0),
        legend.position = c(0.05, 0.95),
        legend.key.width = unit(2.2, "cm"),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = "bold", colour = "black"),
        legend.text.align = 0) +
  guides(color = guide_legend(label.position = "top",
                              label.hjust = 0.5,
                              label.vjust = 1,
                              direction = "horizontal",
                              nrow = 1,
  ))

patch_fig4 <- (trxpt_fpkms/(reg_fpkms| all_juncs | trxptome_juncs)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 52, face = "bold"))

ggsave(plot = patch_fig4, "results/r/figures/figure_4a_d.png",
       width = 40, height = 25, dpi = 350)
