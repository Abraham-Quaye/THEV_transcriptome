# RNA reads coverage and depth of THEV genome
library(tidyverse)
library(glue)
library(plotly)
library(gridExtra)
library(patchwork)
library(ggsci)
library(ggtext)

# load in depth files
# 4hr samples
depth4hrs <- read_tsv("coverage/mappedNCBI_4hrsdepth.txt") %>%
  rename(chr = "#CHROM",
         position = "POS",
         depth = "/Users/abquaye/MyStuff/THEVRNAseqAnalysis/Hisat2_Analyses/mapped/sortedNCBI_4hrsSamples.bam") 


#12hr samples
depth12hrs <- read_tsv("coverage/mappedNCBI_12hrsdepth.txt") %>%
  rename(chr = "#CHROM",
         position = "POS",
         depth = "/Users/abquaye/MyStuff/THEVRNAseqAnalysis/Hisat2_Analyses/mapped/sortedNCBI_12hrsSamples.bam")

# 24hr samples
depth24hrs <- read_tsv("coverage/mappedNCBI_24hrsdepth.txt") %>%
  rename(chr = "#CHROM",
         position = "POS",
         depth = "/Users/abquaye/MyStuff/THEVRNAseqAnalysis/Hisat2_Analyses/mapped/sortedNCBI_24hrsSamples.bam")

# 72 hr samples
depth72hrs <- read_tsv("coverage/mappedNCBI_72hrsdepth.txt") %>% 
  rename(chr = "#CHROM",
         position = "POS",
         depth = "/Users/abquaye/MyStuff/THEVRNAseqAnalysis/Hisat2_Analyses/mapped/sortedNCBI_72hrsSamples.bam")

# plots
(p4 <- depth4hrs %>% 
  ggplot(aes(position,depth)) +
  geom_col(color = "grey") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(1000,26000,1000),
                       labels =  glue("{seq(1,26,1)}Kb")
    ) +
  theme_classic() +
    labs(title = "4hr Timepoint",
         x = element_blank(),
         y = "Coverage Depth") +
    theme(plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
          #plot.margin = margin(rep(20,4)),
          panel.grid.major.y = element_line(linewidth = 0.6, linetype = "dashed"),
          axis.title.y = element_text(size = 16, face = "bold", margin = margin(r=10)),
          axis.text = element_text(size = 12, color = "black")))

(p12 <- depth12hrs %>% 
  ggplot(aes(position,depth)) +
  geom_col(color = "grey50") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(1000,26000,1000),
                       labels =  glue("{seq(1,26,1)}Kb")
    ) +
    labs(title = "12hr Timepoint",
         x = element_blank(),
         y = "Coverage Depth") +
  theme_classic() +
    theme(plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
          #plot.margin = margin(rep(20,4)),
          panel.grid.major.y = element_line(linewidth = 0.6, linetype = "dashed"),
          axis.title.y = element_text(size = 16, face = "bold", margin = margin(r=10)),
          axis.text = element_text(size = 12, color = "black")))

(p24 <- depth24hrs %>% 
        ggplot(aes(position,depth)) +
        geom_col(color = 'grey40') +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(1000,26000,1000),
                       labels =  glue("{seq(1,26,1)}Kb")
    ) +
    labs(title = "24hr Timepoint",
         x = element_blank(),
         y = "Coverage Depth") +
        theme_classic() +
    theme(plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
          #plot.margin = margin(rep(20,4)),
          panel.grid.major.y = element_line(linewidth = 0.6, linetype = "dashed"),
          axis.title.y = element_text(size = 16, face = "bold", margin = margin(r=10)),
          axis.text = element_text(size = 12, color = "black")))

(p72 <- depth72hrs %>% 
        ggplot(aes(position,depth)) +
        geom_col(color = "grey30") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0),
                           breaks = seq(1000,26000,1000),
                           labels =  glue("{seq(1,26,1)}Kb")
                           ) +
    labs(title = "72hr Timepoint",
         x = element_blank(),
         y = "Coverage Depth") +
        theme_classic() +
    theme(plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
          #plot.margin = margin(rep(20,4)),
          panel.grid.major.y = element_line(linewidth = 0.6, linetype = "dashed"),
          axis.title.y = element_text(size = 16, face = "bold", margin = margin(r=10)),
          axis.text = element_text(size = 12, color = "black")))

# multiple plot layouts
# 1. gridExtra package
grid.arrange(p4,p12,p24,p72,
             ncol = 2, nrow = 2)

# 2. patchwork
(p4 / p12 / p24 / p72) + plot_annotation(title = "RNA-seq Mapping Coverage of THEV Genome", tag_levels = "I") &
  theme(plot.tag = element_text(size = 18, face = "bold", hjust = 0, vjust = 0),
        #plot.tag.position = c(1,1),
        plot.title = element_text(face = 'bold', hjust = 0.1, size = 16),
        plot.title.position = "panel")

 
p12 + tableGrob(depth12hrs[16000:16020,])



# try combining all samples
depth4hrs <- depth4hrs%>% 
              mutate(timepoint = "t-4hrs")
depth12hrs <- depth12hrs%>% 
              mutate(timepoint = "t-12hrs")
depth24hrs <- depth24hrs%>% 
              mutate(timepoint = "t-24hrs")
depth72hrs <- depth72hrs%>% 
              mutate(timepoint = "t-72hrs") 

alltimes <- rbind(depth4hrs,depth12hrs,depth24hrs, depth72hrs) %>% 
  mutate(timepoint = factor(timepoint, levels = rev((unique(.$timepoint)))))

alltimes %>% 
  ggplot(aes(position,depth, fill = timepoint, color = timepoint)) +
  geom_col() +
  scale_fill_igv() +
  scale_color_igv() +
  scale_y_log10(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(1000,25000,1000),
                     labels =  glue("{seq(1000,25000,1000)}bp")
  ) +
  labs(title = "All timepoint samples",
       x = element_blank(),
       y = "Coverage Depth") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10))

## coverage data
cov4hrs <- read_tsv("coverage/mappedNCBI_4hrscoverage.txt") %>% 
  mutate(tp = "4hrs") %>% 
  select(numreads, meandepth, tp)

cov12hrs <- read_tsv("coverage/mappedNCBI_12hrscoverage.txt") %>% 
  mutate(tp = "12hrs") %>% 
  select(numreads, meandepth, tp)

cov24hrs <- read_tsv("coverage/mappedNCBI_24hrscoverage.txt") %>% 
  mutate(tp = "24hrs") %>% 
  select(numreads, meandepth, tp)

cov72hrs <- read_tsv("coverage/mappedNCBI_72hrscoverage.txt") %>% 
  mutate(tp = "72hrs") %>% 
  select(numreads, meandepth, tp)

cov_all <- rbind(cov4hrs,cov12hrs,cov24hrs,cov72hrs) %>%
  mutate(tp = factor(tp, levels = c("4hrs", "12hrs", '24hrs', "72hrs")))

lmod <- lm(meandepth ~ numreads, data = cov_all)
slmod <- summary(lmod)


cov_all %>% 
  mutate(tp = recode(tp, `4hrs` = "4hpi", 
                     `12hrs` = "12hpi", 
                     `24hrs` = "24hpi", 
                     `72hrs` = "72hpi")) |> 
  ggplot(aes(numreads,meandepth, color = tp)) +
  geom_smooth(show.legend = F, method = "lm", formula = y ~ x,
              color = 'black', linewidth = 0.2) +
  geom_point(size = 15) +
  geom_textbox(aes(label = glue('**Mean Depth** <br>{round(meandepth,1)}')), 
             nudge_y = 0.25, 
             nudge_x = -0.1,
             size = 6,
             hjust = 0.5,
             show.legend = F, color = "black",
             box.padding = unit(0,"pt"),
             width = NULL,
             box.size = 0) +
  geom_textbox(aes(label = glue('**Reads** <br>{numreads}')), 
             nudge_y = -0.25, 
             nudge_x = 0.2,
             size = 6,
             show.legend = F, color = "black",
             box.padding = unit(0,"pt"),
             width = NULL,
             box.size = 0) +
  annotate(geom = "text", x = 2.5e5, y = 5e2, label = glue("R^2 == {round(slmod$adj.r.squared,4)}"),
           parse = T, size = 12) +
  scale_y_log10() +
  scale_x_log10(limits = c(300,1e8)) +
  labs(title = "Correlation of Mapped Reads to Depth Coverage of THEV genome",
       x = "Total Mapped Reads",
       y = "Mean Coverage Depth",
       color = element_blank()) +
  scale_color_manual(values = rainbow(4)) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.025, color = "black", linetype = "solid"),
        plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b=10)),
        axis.line = element_line(linewidth = 0.4),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size =18, face = "bold"),
        axis.ticks.length.y = unit(0,"cm"),
        legend.text = element_text(size = 12, face = 'bold'), 
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(color = "black", linewidth = 0.1),
        legend.justification = c(0,0),
        legend.position = c(0.05,0.6))
  
# 
# allcovCount <- read.table("THEVmapCoverage.txt",# blank.lines.skip = F, 
#                           skip = 1, comment.char = "", 
#                           header = T,fill = NA) |> 
#  filter(X.rname == "THEV") %>%
#   mutate(numreads = as.numeric(numreads))
# 
# sum(allcovCount$numreads) # total mapped reads from all replicates 
# sum(cov_all$numreads) # total mapped reads from all time-points
#   







