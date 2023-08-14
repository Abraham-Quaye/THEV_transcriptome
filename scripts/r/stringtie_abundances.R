
library(tidyverse)
library(magrittr)
library(ggsci)

# load in stringtie abundance tables
files <- list.files("results/stringtie/", pattern = "t\\d{1,2}S\\d\\.tab",
                    full.names = T) %>% 
  setNames(c("12h.p.i", "12h.p.i", "24h.p.i", "24h.p.i", "24h.p.i", "4h.p.i",
             "4h.p.i", "72h.p.i", "72h.p.i", "72h.p.i"))

stringtie_abun <- map_dfr(files, read_tsv,
                          show_col_types = F,
                          .id = "timepoint") %>%
  dplyr::rename(gene_name = "Gene Name", gene_id = "Gene ID", strand = Strand,
         start = Start, end = End, coverage = Coverage, fpkm = FPKM, tpm = TPM) %>% 
  mutate(gene_name = ifelse(gene_name == "-", "unknown", gene_name)) %>% 
  split(.$timepoint) %>% 
  map(mutate, tot_fpkm_time = sum(fpkm)) %>% 
  do.call("rbind", .) %>% 
  as_tibble() %>%
  select(-c(gene_id)) %>%
  group_by(timepoint, gene_name) %>%
  reframe(count = n(),
          tot_fpkm_time = tot_fpkm_time,
          sum_fpkm = sum(fpkm),
          perc_fpkm = (sum_fpkm / tot_fpkm_time) * 100) %>% 
  unique()

stringtie_abun %>% 
  filter(!gene_name == "unknown") %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("4h.p.i", "12h.p.i", "24h.p.i", "72h.p.i"))) %>% 
  ggplot(aes(gene_name, perc_fpkm, fill = timepoint)) +
  geom_col(position = position_dodge(0.7), width = 0.4) +
  scale_fill_jama() +
  theme_classic()


## plotting values from python script: pile_reads_to_orfs.py


# load in stringtie abundance tables
py_files <- list.files("results/hisat2/bulk", pattern = "subsetTHEV_\\d{1,2}hrs.bam_pileup.csv",
                    full.names = T) %>% 
  setNames(c("12h.p.i", "24h.p.i","4h.p.i", "72h.p.i"))

orf_read_pileups <- map_dfr(py_files, read_csv,
                          show_col_types = F,
                          .id = "timepoint")

orf_read_pileups %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("4h.p.i", "12h.p.i", "24h.p.i", "72h.p.i"))) %>% 
  ggplot(aes(gene_name, log2(count), fill = timepoint)) +
  geom_col(position = position_dodge(0.7), width = 0.4) +
  scale_fill_jama() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.015, 0.015)) +
  labs(title = "Changes in Expression Levels of predicted THEV ORFs Across Four Time Points",
       x = element_blank(),
       y = "Count",
       fill = element_blank()) +
  theme_classic() +
  theme(plot.margin = margin(rep(20, 4)),
        plot.title = element_text(size = 18, colour = "#000000", face = "bold",
                                  hjust = 0.5),
        panel.grid.major.y = element_line(colour = "grey", linewidth = 0.2, linetype = "dashed"),
        panel.grid.minor.y = element_line(colour = "grey", linewidth = 0.2, linetype = "dashed"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "#000000"),
        legend.justification = "center",
        legend.position = "bottom",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.spacing.x = unit(0.5, "cm"))
