#!/usr/bin/env Rscript

library(readxl)
library(ggsci)
library(glue)
library(lubridate)
library(tidyverse)

source("scripts/r/thev_cov_depth.R")

exp2 <- "raw_files/wetlab_data/thev_growthcurve04_2023.xls"
replicate_data <- read_excel(exp2,
                          sheet = "Results",
                          skip = 6,
                          trim_ws = TRUE) %>% 
  dplyr::rename(well = "Well",
                sample_name = "Sample Name",
                conc = "Quantity") %>%
  select(well, sample_name, conc) %>% 
  drop_na(sample_name) %>% 
  mutate(sample_name = str_replace(string = sample_name,
                                   pattern = "(\\d{2})[a-z]{3}(\\w{2})",
                                   replacement = "\\2_\\1")) %>% 
  separate(sample_name, into = c("condition", "hrs_pi"),
           sep = "_") %>% 
  mutate(condition = ifelse(str_detect(condition, "S"), "Inf", "Neg"),
         inf_hrs = paste0(condition, "_", hrs_pi),
         inf_hrs = case_match(inf_hrs,
                              "Inf_00" ~ "Inf_0",
                              "Neg_00" ~ "Neg_0",
                              .default = inf_hrs),
         hrs_pi = as.numeric(hrs_pi)) %>% 
  select(condition, conc, hrs_pi, inf_hrs)


plot_prepped <- replicate_data %>%
  group_by(inf_hrs, hrs_pi) %>% 
  reframe(mean_conc = mean(conc), 
            mean_conc_per_mL = (mean_conc * 10 * 1000),
            replicates = n(),
            stderr_conc_perML = (sd(conc)/sqrt(replicates)) * 10000) %>% 
  mutate(treatment = str_replace(inf_hrs, "([a-zA-Z]{3})_\\d{1,2}", "\\1"))


# add data from first experiment with more time-points (use only extra time-points)
# ---------------------
exp1 <- "raw_files/wetlab_data/thev_growthcurve2021.xls"
rep1_data <- read_excel(exp1,
                        sheet = "Results",
                        skip = 6,
                        trim_ws = TRUE) %>% 
  dplyr::rename(well = "Well",
                sample_name = "Sample Name",
                task = "Task", conc = "Quantity") %>% 
  select(well, task, sample_name, conc) %>% 
  drop_na(sample_name) %>% 
  filter(!str_detect(sample_name, ".*DNA.*")) %>%
  mutate(sample_name = case_when(well == "B3" ~ "Inf1-D0",
                                 well == "B4" ~ "Inf2-D0",
                                 TRUE ~ sample_name)) %>%
  mutate(sample_name = str_replace(sample_name,
                                   "([a-zA-Z]+)\\d?-D(\\d)",
                                   "\\1_\\2")) %>% 
  separate_wider_delim(sample_name,
                       names = c("condition", "day"),
                       delim = "_") %>%
  mutate(inf_day = paste0(condition, "_", day),
         day = as.numeric(day),
         hrs_pi = day * 24,
         inf_hrs = paste0(condition, "_", hrs_pi)) %>% 
  select(condition, conc, hrs_pi, inf_hrs) %>% 
  filter(hrs_pi > 72)

fulltimes <- bind_rows(replicate_data, rep1_data)


full_prepped <- fulltimes %>%
group_by(condition, hrs_pi) %>% 
  reframe(mean_conc = mean(conc), 
            mean_conc_per_mL = (mean_conc * 10 * 1000),
            replicates = n(),
            stderr_conc_perML = (sd(conc)/sqrt(replicates)) * 10000)

growth_curve <- full_prepped %>% 
  ggplot(aes(hrs_pi, mean_conc_per_mL, color = condition)) +
  geom_errorbar(aes(ymax = mean_conc_per_mL + stderr_conc_perML,
                    ymin = mean_conc_per_mL - stderr_conc_perML),
                width = 0.8, show.legend = F,
                linewidth = 1, color = "#000000") +
  geom_point(size = 5) +
  geom_line(linewidth  = 1.5) +
  scale_color_aaas(breaks = c("Inf", "Neg"),
                   labels = c("Infected", "Mock-Infected")) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     breaks = full_prepped$hrs_pi,
                     labels = c(paste0(full_prepped$hrs_pi, "hpi"))) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(title = "THEV Growth Curve in RP-19 Turkey B-Cells",
       x = element_blank(),
       y = "Virus Titer (GCN/mL)",
       size = element_blank()
  ) +
  theme(plot.background = element_rect(fill = "#ffffff"),
        plot.margin = margin(rep(5,4)),
        panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.y = element_line(linewidth = 0.2,
                                          color = "grey",
                                          linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 28, face = "bold", colour = "#000000",
                                  hjust = 0.5, margin = margin(b = 10)),
        plot.caption = element_text(colour = "#000000", size = 12,
                                    hjust = 0, face = "italic"),
        plot.caption.position = "plot",
        axis.line = element_line(linewidth = 0.4, colour = "#000000"),
        axis.text.x = element_text(size = 18, color = "#000000", face = "bold"),
        axis.text.y = element_text(size = 14, color = "#000000", face = "bold"),
        axis.title.y = element_text(size = 22, face = "bold", colour = "#000000",
                                    margin = margin(r = 10)),
        axis.ticks.length.y = unit(0, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, face = "bold", colour = "#000000",
                                   margin = margin(rep(5, 4))),
        legend.key = element_rect(fill = "#ffffff",
                                  linewidth = 5),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(1,"cm"),
        legend.background = element_rect(fill = "#ffffff",
                                         color = "#000000",
                                         linewidth = 0.1),
        legend.justification = c(0, 0),
        legend.position = c(0.1, 0.7)
  )


fig_2 <- (comp_all/ growth_curve) +
  plot_layout(heights = c(1.4, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 26, face = "bold"))


ggsave(plot = fig_2,filename = "results/r/figures/fig_2.png",
       width = 14, height = 16, dpi = 500)

