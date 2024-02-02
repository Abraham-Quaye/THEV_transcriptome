library(tidyverse)
library(magrittr)
library(ggsci)

# load files
orf_reads_files <- list.files("results/hisat2/bulk",
           pattern = "subsetTHEV_\\d{1,2}hrs\\.bam_pileup.csv",
           full.names = T) %>%
  magrittr::set_names(c("12hpi", "24hpi", "4hpi", "72hpi"))


orf_read_counts <- map_dfr(orf_reads_files, read.csv, header = T, .id = "timepoint") %>%
  mutate(timepoint = factor(timepoint, levels = c("4hpi", "12hpi", "24hpi", "72hpi")),
         gene_name = factor(gene_name, levels = c("ORF1", "Hyd", "IVa2", "AdPOL",
                                                  "pTP", "DBP", "UXP", "100K", "22K",
                                                  "33K", "pVIII", "E3", "ORF8", "52K",
                                                  "pIIIa", "III", "pVII", "pX", "pVI",
                                                  "Hexon", "Protease", "Fiber", "ORF7")))

orf_read_counts %>%
  filter(timepoint != "4hpi") %>% 
  ggplot(aes(gene_name, tpm, group = timepoint, fill = timepoint)) +
  geom_col(position = "dodge") +
  scale_fill_jco() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Read Coverage of the Predicted ORFs of THEV",
       x = "Predicted ORF",
       y = "TPM") +
  theme_classic() +
  theme(plot.background = element_rect(fill = "#ffffff"),
        plot.margin = margin(rep(30,4)),
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
        axis.text.x = element_text(size = 15, color = "#000000", face = "bold"),
        axis.text.y = element_text(size = 14, color = "#000000", face = "bold"),
        axis.title = element_text(size = 22, face = "bold", colour = "#000000",
                                    margin = margin(r = 10, t = 10)),
        axis.ticks.length.y = unit(0, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, face = "bold", colour = "#000000",
                                   margin = margin(rep(5, 4))),
        legend.key = element_rect(fill = "#ffffff", colour = "#ffffff"),
        legend.key.width = unit(1,"cm"),
        legend.key.height = unit(1,"cm"),
        legend.background = element_rect(fill = "#ffffff"),
        legend.justification = c(0, 0),
        legend.position = c(0.05, 0.75)
  )
