#!/usr/bin/env Rscript


library(magrittr)
library(glue)
library(patchwork)
library(ggsci)
library(ggtext)
library(grDevices)
library(tidyverse)

source("scripts/r/thev_genomic_map.R")

### remove older files before generating new files
path <- "results/r/figures"
files <- list.files(path, pattern = ".+(.pdf|.png|.jpg|.jpeg)", 
                    full.names = TRUE, recursive = TRUE)
unlink(files, recursive = FALSE, force = FALSE)

# load in depth files

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
  
# -------------------
## split the table into individual time-points and plot iteratively with map()
each_plot <- all_depths %>% 
  split(.$timepoint) %>% 
  map(~ggplot(., aes(position, depth)) +
        geom_col(color = glue("{.$color}")) +
        labs(title = glue("{.$titles}"),
             x = element_blank(),
             y = "Mapping Depth") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0),
                           breaks = seq(1000, 26000, 2000),
                           labels = glue("{seq(1,26,2)}kb")) +
        theme_classic() +
        theme(plot.title = element_text(size = 18,
                                        face = "bold",
                                        hjust = 0.5,
                                        margin = margin(t = 10)),
              panel.grid.major.y = element_line(linewidth = 0.6,
                                                linetype = "dashed"),
              axis.title.y = element_text(size = 14,
                                          face = "bold",
                                          margin = margin(r = 10, l = 10)),
              axis.text.y = element_text(size = 10, color = "black"),
              axis.text.x = element_text(size = 8.5, color = "black", face = "bold")
            ))

# save plot for each time-point

# ggsave("depth_4hrs.png",
#        plot = each_plot$`4hpi`, path = "results/r/figures",
#        width = 10, height = 5)
# 
# ggsave("depth_12hrs.png",
#        plot = each_plot$`12hpi`, path = "results/r/figures",
#        width = 10, height = 5)
# 
# ggsave("depth_24hrs.png",
#        plot = each_plot$`24hpi`, path = "results/r/figures",
#        width = 10, height = 5)
# 
# ggsave("depth_72hrs.png",
#        plot = each_plot$`72hpi`, path = "results/r/figures",
#        width = 10, height = 5)

# genomic map plot
genome <- make_genomic_map("raw_files/annotations/THEVannotated_genesOnly.bed")


# 2. patchwork
p_alltime <- (each_plot$`4hpi`/each_plot$`12hpi`/each_plot$`24hpi`/each_plot$`72hpi`/genome) +
  plot_annotation(title = "RNA-seq Mapping Depth of THEV Genome") +
  plot_layout(heights = c(rep(6, 4), 2)) &
  theme(plot.tag = element_text(
    size = 18,
    face = "bold",
    hjust = 0,
    vjust = 0),
    plot.title = element_text(face = "bold",
                              hjust = 0.1,
                              size = 16),
    plot.title.position = "panel")

ggsave("patch_alltimes.png",
       plot = p_alltime, path = "results/r/figures",
       width = 12, height = 14, dpi = 1000)

## --------------------

comp_all <- all_depths %>%
  select(timepoint, position, depth) %>% 
  mutate(depth = ifelse(depth < 1, 1, depth)) %>% 
  ggplot(aes(position, sqrt(depth), group = timepoint, fill = timepoint, color = timepoint)) +
  geom_area(key_glyph = "crossbar") +
  scale_fill_igv(guide = guide_legend(keywidth = unit(1, "cm"))) +
  scale_color_igv() +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(1000, 26000, 1000),
                     labels = glue("{seq(1,26,1)}kb")) +
  coord_trans(clip = "off", expand = F) +
  labs(title = "RNA-seq Mapping Depth of THEV Genome",
       x = element_blank(),
       y = expression(sqrt("Mapping Depth")),
       fill = element_blank(),
       color = element_blank()) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10,
                                   face = "bold",
                                   color = "#000000",
                                   margin = margin(l = 10, r = 5)),
        axis.text.x = element_text(size = 10,
                                   face = "bold",
                                   color = "#000000",
                                   margin = margin(b = 10, t = 5)),
        plot.title = element_text(size = 20, 
                                  face = "bold", 
                                  hjust = 0.5,
                                  margin = margin(t = 10)),
        axis.title.y = element_text(size = 22, face = "bold",
                                    margin = margin(l = 10)),
        legend.justification = c(0, 0),
        legend.position = c(0.04, 0.7),
        legend.margin = margin(rep(10, 4)),
        legend.background = element_rect(color = "grey", linewidth = 0.2),
        legend.text = element_text(size = 14, 
                                   margin = margin(rep(10, 4)),
                                   face = "bold"))


compare_all <- (comp_all/genome) +
  plot_layout(heights = c(11, 1)) 

ggsave("overlay_alltimes.png",
       plot = compare_all, path = "results/r/figures",
       width = 15, height = 10, dpi = 1000)


################# coverage data #########################

### find all the coverage files
all_cov_files <- list.files("results/hisat2/coverage",
                            pattern = "host_thev_cov\\d+.txt",
                            full.names = TRUE) %>% 
  setNames(c("12hpi", "24hpi", "4hpi", "72hpi"))

# read depth_files as one master tibble
all_covs <- map_dfr(all_cov_files, read_tsv,
                    show_col_types = FALSE,
                    comment = "Coverage",
                    .id = "timepoint") %>%
  rename("organism" = "#rname") %>% 
  mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi"))) %>% 
  map_at(c(3:10), as.numeric) %>% 
  as_tibble() %>%
  mutate(organism = ifelse(organism == "AY849321.1", "thev", "m.gallopavo")) %>% 
  group_by(organism, timepoint) %>%
  reframe(total_mapped = sum(numreads),
          mean_depth = mean(meandepth), mean_cov = mean(coverage))

cov_thev <- all_covs %>% filter(organism == "thev")


lmod <- lm(mean_depth ~ total_mapped, 
           data = cov_thev %>% select(total_mapped, mean_depth))
slmod <- summary(lmod)


corr <- cov_thev %>%
  ggplot(aes(total_mapped, mean_depth, color = timepoint)) +
  geom_smooth(show.legend = FALSE, se = FALSE,
              method = "lm", 
              formula = y ~ x,
              color = "black",
              linewidth = 0.2) +
  geom_point(size = 15) +
  geom_text(aes(label = glue("Mean Depth \n{round(mean_depth,1)}")),
            nudge_y = 0.25,
            nudge_x = -0.05,
            size = 4,
            fontface = "bold",
            show.legend = F, color = "black") +
  geom_text(aes(label = glue("Reads \n{total_mapped}")),
            fontface = "bold",
            nudge_y = -0.25,
            nudge_x = 0.1,
            size = 4,
            show.legend = F, color = "black") +
  annotate(geom = "text", x = 2.5e5, y = 5e2, 
           label = glue("R^2 == {round(slmod$adj.r.squared,4)}"),
           parse = T, size = 12) +
  scale_y_log10() +
  scale_x_log10(limits = c(300, 1e8)) +
  labs(title = "Correlation of Mapped Reads to Coverage Depth of THEV genome",
       x = "Total Mapped Reads",
       y = "Mean Depth/Basepair",
       color = element_blank()) +
  scale_color_manual(values = rainbow(4)) +
  theme(plot.background = element_blank(),
        plot.margin = margin(rep(20,4)),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.01, 
                                          color = "grey80", 
                                          linetype = "dashed"),
        plot.title = element_text(size = 28, face = "bold", 
                                  hjust = 0.5, margin = margin(b = 10)),
        axis.line = element_line(linewidth = 0.4),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 22, face = "bold"),
        axis.ticks.length.y = unit(0, "cm"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(color = "black", linewidth = 0.1),
        legend.justification = c(0, 0),
        legend.position = c(0.05, 0.6))

ggsave("correlate_alltimes.png",
       plot = corr, path = "results/r/figures",
       width = 20, height = 14, dpi = 1000)