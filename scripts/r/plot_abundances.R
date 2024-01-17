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
          axis.text.x = element_text(size = 28, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"),
          legend.justification = c(0, 0),
          legend.position = c(0.05, 0.9),
          legend.key.width = unit(2, "cm"),
          legend.title = element_blank(),
          legend.text = element_text(size = 15, face = "bold", colour = "black"),
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

# patch_expr <- (all_juncs | trxptome_juncs)

# ggsave("junc_abundances.png",
#        plot = patch_expr, path = "results/r/figures",
#        width = 20, height = 12, dpi = 500)


# -------------------
## splice donor and acceptor frequencies

# acceptors and donors in trxptome
bulk_trxptome_ss_seq <- count_trxptome_ss(unq_bulk_juncs) %>%
  distinct(junc_ss, junc_ts, .keep_all = T) %>%
  dplyr::group_by(splice_site) %>%
  reframe(freq = n()) %>% 
  mutate(sum_junc_count = sum(freq),
         percent_abund = (freq / sum_junc_count) * 100)

# plot trxptome acceptors and donor frequencies
ss <- bulk_trxptome_ss_seq %>%
  ggplot(aes(splice_site, percent_abund, color = splice_site)) +
  geom_point(show.legend = F, size = 10) +
  geom_segment(aes(x = splice_site, xend = splice_site,
                   y = 0, yend = percent_abund),
               linewidth = 2,
               show.legend = F) +
  geom_text(aes(label = glue("{round(percent_abund, 1)}%")),
             nudge_y = 3.5,
            size = 10, fontface = "bold", color = "#000000") +
  labs(x = "Splice Site Donor-Acceptor",
       y = "Frequency") +
  scale_color_lancet() +
  scale_y_continuous(expand = c(0.01,0.01),
                     labels = scales::label_percent(scale = 1)) +
  scale_x_discrete(expand = c(0.13, 0.13)) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(plot.margin = margin(rep(30, 4)),
        panel.grid.major.y = element_line(linewidth = 0.4, color = "grey",
                                          linetype = "dashed"),
        axis.title = element_text(size = 16, face = "bold", margin = margin(r = 15, t = 25)),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(color = "grey40"),
        axis.ticks.length = unit(0, "pt")
        )

# All acceptors and donors
# all_ss_seq <- bulk_junc_stats %>%
#   distinct(timepoint, start, end, .keep_all = T) %>%
#   mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>% 
#   group_by(timepoint, splice_site) %>%
#   reframe(tally = n(),
#           total_reads = sum(read_count)) %>%
#   split(.$timepoint) %>%
#   map(mutate, tot_tally_tp = sum(tally)) %>% 
#   do.call("rbind", .) %>%
#   mutate(percent_abund = (tally / tot_tally_tp) * 100) %>% 
#   ggplot(aes(splice_site, percent_abund, group = timepoint, color = timepoint)) +
#   geom_point(show.legend = F, size = 10) +
#   geom_segment(aes(x = splice_site, xend = splice_site,
#                    y = 0, yend = percent_abund),
#                linewidth = 2,
#                show.legend = F) +
#   geom_text(aes(label = glue("{round(percent_abund, 1)}%")),
#             nudge_y = 4.5,
#             size = 5, fontface = "bold", color = "#000000") +
#   facet_wrap(~ timepoint, scales = "free") +
#   labs(x = "Splice Site Donor-Acceptor",
#        y = "Frequency") +
#   scale_color_lancet() +
#   scale_y_continuous(expand = c(0.01,0.01),
#                      limits = c(0, 105),
#                      labels = scales::label_percent(scale = 1)) +
#   scale_x_discrete(expand = c(0.05, 0.05)) +
#   coord_cartesian(clip = "off") +
#   theme_classic() +
#   theme(plot.margin = margin(rep(30, 4)),
#         panel.grid.major.y = element_line(linewidth = 0.4, color = "grey",
#                                           linetype = "dashed"),
#         axis.title = element_text(size = 16, face = "bold", margin = margin(r = 15, t = 15)),
#         axis.text.x = element_text(size = 16, face = "bold"),
#         axis.text.y = element_text(size = 14, face = "bold"),
#         axis.line = element_line(color = "grey40"),
#         axis.ticks.length = unit(0, "pt"),
#         strip.background = element_rect(linewidth = 0.1, fill = "grey"),
#         strip.text.x = element_text(size = 12, face = "bold", margin = margin(0.2,0,0.2,0, "cm")),
#         strip.placement = "outside",
#         panel.spacing = unit(1, "cm")
#   )

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
        axis.text.x = element_text(size = 28, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, face = "bold"),
        legend.justification = c(0, 0),
        legend.position = c(0.02, 0.8),
        legend.key.size = unit(1.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, face = "bold", colour = "black"),
        legend.text.align = 0)

# ggsave(plot = trxpt_fpkms, "results/r/figures/trxpt_fpkm_percent_abund.png",
#        width = 18, height = 12, dpi = 500)


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
        legend.position = c(0.05, 0.9),
        legend.key.width = unit(2, "cm"),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, face = "bold", colour = "black"),
        legend.text.align = 0) +
  guides(color = guide_legend(label.position = "top",
                              label.hjust = 0.5,
                              label.vjust = 1,
                              direction = "horizontal",
                              nrow = 1,
  ))

# ggsave(plot = reg_fpkms, "results/r/figures/region_fpkm_percent_abund.png",
#        width = 18, height = 12, dpi = 500)

# plot for figure 4A-D


patch_fig4 <- (trxpt_fpkms/(reg_fpkms| all_juncs | trxptome_juncs)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 52, face = "bold"))

ggsave(plot = patch_fig4, "results/r/figures/figure_4a_d.png",
       width = 45, height = 30, dpi = 350)
# ------------------------------
# distribution of fpkm levels in samples

# fpkm_dist <- t_exp_levels %>%
#   mutate(fpkm = log2(fpkm + 1),
#          timepoint = factor(timepoint,
#                             levels = c("4h.p.i", "12h.p.i", "24h.p.i", "72h.p.i"))) %>%
#   ggplot(aes(timepoint, fpkm, fill = timepoint)) +
#   geom_boxplot() +
#   geom_jitter(show.legend = F, width = 0.25, size = 3, alpha = 0.6) +
#   scale_y_continuous(expand = c(0.01, 0.01)) +
#   scale_fill_igv() +
#   labs(title = "Distribution of Transcript FPKM Values Across Four Timepoints",
#        x = "Time-Point",
#        y = "log2(FPKM + 1)",
#        fill = element_blank()) +
#   theme_classic() +
#   theme(plot.margin = margin(rep(20, 4)),
#         plot.title = element_text(size = 18, colour = "#000000", face = "bold",
#                                   hjust = 0.5),
#         axis.title = element_text(size = 16, face = "bold"),
#         axis.text = element_text(size = 14, color = "#000000"),
#         legend.justification = c(0, 0),
#         legend.position = c(0.92, 0.85),
#         legend.key.size = unit(1, "cm"),
#         legend.text = element_text(size = 14, face = "bold"))
# 
# ggsave(plot = fpkm_dist, "results/r/figures/fpkm_dist_by_time.png",
#        width = 14, height = 10, dpi = 500)
