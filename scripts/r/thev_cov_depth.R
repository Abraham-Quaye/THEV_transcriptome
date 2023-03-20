#!/usr/bin/env Rscript


library(tidyverse)
library(glue)
library(gridExtra)
library(patchwork)
library(ggsci)
library(ggtext)

### remove older files before generating new files
path <- "results/r/figures"
files <- list.files(path, pattern = ".+(.pdf|.png|.jpg|.jpeg)", 
                    full.names = TRUE, recursive = TRUE)
unlink(files, recursive = FALSE, force = FALSE)

# load in depth files

### find all the depth files and read them as one master dataframe/tibble
all_depth_files <- list.files("results/hisat2/coverage",
                              pattern = "[a-z_]+\\d+[a-z]+\\.txt",
                              full.names = TRUE)
names(all_depth_files) <- c("12hpi", "24hpi", "4hpi", "72hpi")

all_depths <- map_dfr(all_depth_files, read_tsv,
                      col_names = FALSE,
                      comment = "#",
                      .id = "timepoint") %>% 
  mutate(timepoint = factor(timepoint,
                            levels = c("72hpi",
                                       "24hpi",
                                       "12hpi",
                                       "4hpi")))

colnames(all_depths) <- c("timepoint","genome", "position", "depth")
  

## split the table into individual time-points to plot separately
each_timepoint <- all_depths %>% 
  split(.$timepoint)

# plots
looks <- (theme_classic() +
            theme(plot.title = element_text(size = 16,
                                            face = "bold",
                                            hjust = 0.5), 
                  panel.grid.major.y = element_line(linewidth = 0.6,
                                                    linetype = "dashed"),
                  axis.title.y = element_text(size = 16,
                                              face = "bold",
                                              margin = margin(r = 10)),
                  axis.text.y = element_text(size = 10, color = "black"),
                  axis.text.x = element_text(size = 8.5, color = "black", face = "bold")
            ))

# plot 4hpi mapping depth

p4 <- each_timepoint$`4hpi` %>% 
  ggplot(aes(position, depth)) +
  geom_col(color = "grey") +
  labs(title = "4hrs Post-infection",
       x = element_blank(),
       y = "Mapping Depth") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(1000, 26000, 1000),
                     labels = glue("{seq(1,26,1)}kb")) +
  looks

ggsave("depth_4hrs.pdf",
       plot = p4, path = "results/r/figures",
       width = 10, height = 5)


# plot 12hpi mapping depth

p12 <- each_timepoint$`12hpi` %>% 
  ggplot(aes(position, depth)) +
  geom_col(color = "grey50") +
  labs(title = "12hrs Post-infection",
       x = element_blank(),
       y = "Mapping Depth") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(1000, 26000, 1000),
                     labels = glue("{seq(1,26,1)}kb")) +
  looks

ggsave("depth_12hrs.pdf",
       plot = p12, path = "results/r/figures",
       width = 10, height = 5)


# plot 24hpi mapping depth

p24 <- each_timepoint$`24hpi` %>% 
  ggplot(aes(position, depth)) +
  geom_col(color = "grey40") +
  labs(title = "24hrs Post-infection",
       x = element_blank(),
       y = "Mapping Depth") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(1000, 26000, 1000),
                     labels = glue("{seq(1,26,1)}kb")) +
  looks

ggsave("depth_24hrs.pdf",
       plot = p24, path = "results/r/figures",
       width = 10, height = 5)


# plot 72hpi mapping depth

p72 <- each_timepoint$`72hpi` %>% 
  ggplot(aes(position, depth)) +
  geom_col(color = "grey30") +
  labs(title = "72hrs Post-infection",
       x = element_blank(),
       y = "Mapping Depth") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(1000, 26000, 1000),
                     labels = glue("{seq(1,26,1)}kb")) +
  looks

ggsave("depth_72hrs.pdf",
       plot = p72, path = "results/r/figures",
       width = 10, height = 5)


# multiple plot layouts
# 1. gridExtra package
# alltime <- grid.arrange(p4, p12, p24, p72,
#             ncol = 1, nrow = 4)
# ggsave("grid_alltimes.pdf",
#   plot = alltime, path = "results/r/figures",
#   width = 5, height = 10
# )


# 2. patchwork
p_alltime <- (p4 / p12 / p24 / p72) +
  plot_annotation(title = "RNA-seq Mapping Coverage of THEV Genome",
                  tag_levels = "I") &
  theme(plot.tag = element_text(
    size = 18,
    face = "bold",
    hjust = 0,
    vjust = 0),
    plot.title = element_text(face = "bold",
                              hjust = 0.1,
                              size = 16),
    plot.title.position = "panel")

ggsave("patch_alltimes.pdf",
       plot = p_alltime, path = "results/r/figures",
       width = 10, height = 15)



compare_all <- all_depths %>%
  ggplot(aes(position, depth, fill = timepoint, color = timepoint)) +
  geom_col() +
  scale_fill_igv(guide = guide_legend(
    keywidth = unit(1, "cm"),
    keyheight = unit(1.2, "cm"))) +
  scale_color_igv() +
  scale_y_sqrt(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(1000, 26000, 1000),
                     labels = glue("{seq(1,26,1)}kb")) +
  labs(title = "All timepoint samples",
       x = element_blank(),
       y = "Mapping Depth",
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
        axis.title.y = element_text(size = 18, face = "bold",
                                    margin = margin(l = 10)),
        legend.justification = c(0, 0),
        legend.position = c(0.04, 0.75),
        legend.margin = margin(rep(10, 4)),
        legend.background = element_rect(color = "grey", linewidth = 0.2),
        legend.text = element_text(size = 14, 
                                   margin = margin(rep(10, 4)),
                                   face = "bold"))

ggsave("overlay_alltimes.pdf",
       plot = compare_all, path = "results/r/figures",
       width = 15, height = 10)


# coverage data


cov_all <- read_tsv("results/hisat2/coverage/bulk_coverage.txt",
                    comment = "Coverage", show_col_types = FALSE) %>% 
  rename("organism" = "#rname") %>% 
  filter(organism != "#rname") %>% 
  mutate(timepoint = c("4hpi", "24hpi", "12hpi","72hpi"),
         timepoint = factor(timepoint,levels = c("4hpi","12hpi","24hpi","72hpi"))) %>% 
  map_at(c(2:9), as.numeric) %>% 
  as_tibble()


lmod <- lm(meandepth ~ numreads, 
           data = cov_all %>% select(numreads, meandepth))
slmod <- summary(lmod)


corr <- cov_all %>%
  ggplot(aes(numreads, meandepth, color = timepoint)) +
  geom_smooth(show.legend = FALSE, se = FALSE,
              method = "lm", 
              formula = y ~ x,
              color = "black",
              linewidth = 0.2) +
  geom_point(size = 15) +
  geom_text(aes(label = glue("Mean Depth \n{round(meandepth,1)}")),
            nudge_y = 0.25,
            nudge_x = -0.1,
            size = 6,
            fontface = "bold",
            show.legend = F, color = "black") +
  geom_text(aes(label = glue("Reads \n{numreads}")),
            fontface = "bold",
            nudge_y = -0.25,
            nudge_x = 0.2,
            size = 6,
            show.legend = F, color = "black") +
  annotate(geom = "text", x = 2.5e5, y = 5e2, 
           label = glue("R^2 == {round(slmod$adj.r.squared,4)}"),
           parse = T, size = 12) +
  scale_y_log10() +
  scale_x_log10(limits = c(300, 1e8)) +
  labs(title = "Correlation of Mapped Reads to Coverage Depth of THEV genome",
       x = "Total Mapped Reads",
       y = "Mean Depth",
       color = element_blank()) +
  scale_color_manual(values = rainbow(4)) +
  theme(plot.background = element_blank(),
        plot.margin = margin(rep(20,4)),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.01, 
                                          color = "grey", 
                                          linetype = "solid"),
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

ggsave("correlate_alltimes.pdf",
       plot = corr, path = "results/r/figures",
       width = 20, height = 14)
