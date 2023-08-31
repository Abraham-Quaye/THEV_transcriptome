#!/usr/bin/env Rscript

library(ballgown)
library(devtools)
library(genefilter)
library(ggsci)
library(tidyverse)

# make experimental data dataframe
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
  pivot_longer(cols = starts_with("FPKM"),
               names_to = "samples",
               values_to = "fpkm") %>%
  mutate(timepoint = case_when(str_detect(samples, "_4hrs") ~ "4h.p.i",
                               str_detect(samples, "_12hrs") ~ "12h.p.i",
                               str_detect(samples, "_24hrs") ~ "24h.p.i",
                               str_detect(samples, "_72hrs") ~ "72h.p.i",
                               TRUE ~ NA_character_),
         timepoint = factor(timepoint,
                            levels = c("4h.p.i", "12h.p.i", "24h.p.i", "72h.p.i"))) %>%
  select(c("timepoint", "t_name", region = "gene_name", "strand", "start", "end", "fpkm")) %>%
  split(.$timepoint) %>%
  map(mutate, tot_fpkm_tp = sum(fpkm)) %>% 
  do.call("rbind", .)


# ----------------
# expression levels by region
t_exp_lev_byregion <- t_exp_levels %>%
  group_by(timepoint, region, .drop = F) %>%
  reframe(fpkm_reg = sum(fpkm),
          tot_fpkm_tp = mean(tot_fpkm_tp),
          fpkm_reg_percent = round((fpkm_reg / tot_fpkm_tp) * 100, 2)) %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("4h.p.i", "12h.p.i", "24h.p.i", "72h.p.i")))


# plot transcript expression levels
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
        axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 15)),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        legend.justification = c(0, 0),
        legend.position = c(0.05, 0.9),
        legend.key.width = unit(2, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, face = "bold", colour = "black"),
        legend.text.align = 0) +
  guides(color = guide_legend(label.position = "top",
                              label.hjust = 0.5,
                              label.vjust = 1,
                              direction = "horizontal",
                              nrow = 1,
  ))

ggsave(plot = reg_fpkms, "results/r/figures/region_fpkm_percent_abund.png",
       width = 18, height = 12, dpi = 500)
## ---------------------
# expression levels by transcript
t_expr_each <- t_exp_levels %>%
  group_by(timepoint, t_name) %>%
  reframe(fpkm_trxpt = sum(fpkm),
          tot_fpkm_tp = mean(tot_fpkm_tp),
          fpkm_trxpt_percent = round((fpkm_trxpt / tot_fpkm_tp) * 100, 2)) %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("4h.p.i", "12h.p.i", "24h.p.i", "72h.p.i")),
         t_name = dplyr::case_match(t_name,
                                    "trxptA_novel" ~ "E1: Hyd_iso_1",
                                    "trxptB_ORF1" ~ "E1: ORF1_novel_iso",
                                    "trxptC_hyd" ~ "E1: Hyd_iso_2",
                                    "trxptD_ORF4" ~ "E1: ORF4_novel",
                                    "DBP" ~ "E2A: DBP",
                                    "ptp_pol1" ~ "E2B: pTP/Pol_iso_1",
                                    "ptp_pol2" ~ "E2B: pTP/Pol_iso_2",
                                    "trunc_hypothetical" ~ "E2B: hypo_trunc",
                                    "IVa2" ~ "IM: IVa2",
                                    "100K:Fib" ~ "E3: 100K-Fiber",
                                    "33K_pVIII_E3_Fib" ~ "E3: 33K-Fiber",
                                    "33K_pVIII_E3_lngr_tts" ~ "E3: 33K-E3_iso_1",
                                    "33K_pVIII_E3_xtra_exon" ~ "E3: 33K-E3_iso_2",
                                    "22K:Fib" ~ "E3: 22K-Fiber",
                                    "E3_Fib_ORF7" ~ "E3: E3-ORF7",
                                    "ORF8" ~ "E4: ORF8",
                                    "33K:E3" ~ "MLP: 33K-E3",
                                    "52K" ~ "MLP: 52K",
                                    "Hexon" ~ "MLP: Hexon",
                                    "Hexon:protease" ~ "MLP: Hexon/Protease",
                                    "TPL_exons_trunc" ~ "MLP: truncated_TPL",
                                    "MLP_Fib_ORF7" ~ "MLP: Fiber/ORF7",
                                    "pIIIa:pVII" ~ "MLP: pIIIa-pVII",
                                    "pVII:protease" ~ "MLP: pVII-Protease_iso_1",
                                    "pVII:protease_2" ~ "MLP: pVII-Protease_iso_2",
                                    "III_pVII" ~ "MLP: III/pVII",
                                    "pVII" ~ "MLP: pVII",
                                    "pX:Hexon" ~ "MLP: pX-Hexon",
                                    .default =  t_name
         )) %>%
  mutate(t_name = factor(t_name,
                         levels = c("E1: Hyd_iso_1", "E1: ORF1_novel_iso",
                                    "E1: Hyd_iso_2", "E1: ORF4_novel",
                                    "E2A: DBP",  "E2B: pTP/Pol_iso_1",
                                    "E2B: pTP/Pol_iso_2", "E2B: hypo_trunc",
                                    "IM: IVa2", "E3: 100K-Fiber",
                                    "E3: 33K-Fiber", "E3: 33K-E3_iso_1",
                                    "E3: 33K-E3_iso_2", "E3: 22K-Fiber",
                                    "E3: E3-ORF7", "E4: ORF8", "MLP: 33K-E3",
                                    "MLP: 52K", "MLP: Hexon", "MLP: Hexon/Protease",
                                    "MLP: truncated_TPL", "MLP: Fiber/ORF7",
                                    "MLP: pIIIa-pVII", "MLP: pVII-Protease_iso_1",
                                    "MLP: pVII-Protease_iso_2", "MLP: III/pVII",
                                    "MLP: pVII", "MLP: pX-Hexon")))
# plot
trxpt_fpkms <- t_expr_each %>%
  filter(timepoint != "4h.p.i") %>%
  ggplot(aes(t_name, fpkm_trxpt_percent, group = timepoint, fill = timepoint)) +
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
        axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 15)),
        axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14, face = "bold"),
        legend.justification = c(0, 0),
        legend.position = c(0.9, 0.8),
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, face = "bold", colour = "black"),
        legend.text.align = 0)

ggsave(plot = trxpt_fpkms, "results/r/figures/trxpt_fpkm_percent_abund.png",
       width = 18, height = 12, dpi = 500)

# ------------------------------
# distribution of fpkm levels in samples

fpkm_dist <- t_exp_levels %>%
  mutate(fpkm = log2(fpkm + 1),
         timepoint = factor(timepoint,
                            levels = c("4h.p.i", "12h.p.i", "24h.p.i", "72h.p.i"))) %>%
  ggplot(aes(timepoint, fpkm, fill = timepoint)) +
  geom_boxplot() +
  geom_jitter(show.legend = F, width = 0.25, size = 3, alpha = 0.6) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  scale_fill_igv() +
  labs(title = "Distribution of Transcript FPKM Values Across Four Timepoints",
       x = "Time-Point",
       y = "log2(FPKM + 1)",
       fill = element_blank()) +
  theme_classic() +
  theme(plot.margin = margin(rep(20, 4)),
        plot.title = element_text(size = 18, colour = "#000000", face = "bold",
                                  hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "#000000"),
        legend.justification = c(0, 0),
        legend.position = c(0.92, 0.85),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 14, face = "bold"))

ggsave(plot = fpkm_dist, "results/r/figures/fpkm_dist_by_time.png",
       width = 14, height = 10, dpi = 500)





 
