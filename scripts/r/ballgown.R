

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

# load in annotated ORFs from NCBI
pred_orfs <- gffReadGR("raw_files/annotations/thev_predicted_genes.gtf",
                       splitByTranscript = T)

# annotate assembled transcripts with predicted ORFs from NCBI
assembled_trxpts <- structure(raw_data)$trans
trxpt_anno <- annotate_assembly(assembled_trxpts, pred_orfs)

ballgown::plotMeans("XLOC_000001", raw_data, groupvar = "timepoint",
                    groupname = "all" ,labelTranscripts = T)


structure(raw_data)$exon

ex_exp_levels <- eexpr(raw_data)

# plot transcript expression levels
t_exp_levels <- texpr(raw_data, meas = "all")

t_exp_levels %>%
  as_tibble() %>%
  pivot_longer(cols = starts_with("FPKM"),
               names_to = "samples",
               values_to = "fpkm") %>%
  mutate(fpkm = log2(fpkm + 1),
         samples = factor(samples,
                          levels = c("FPKM.abund_4hrsS1", "FPKM.abund_4hrsS2",
                                     "FPKM.abund_12hrsS1", "FPKM.abund_12hrsS3",
                                     "FPKM.abund_24hrsS1", "FPKM.abund_24hrsS2",
                                     "FPKM.abund_24hrsS3", "FPKM.abund_72hrsS1",
                                     "FPKM.abund_72hrsS2", "FPKM.abund_72hrsS3")),
         timepoint = case_when(str_detect(samples, "_4hrs") ~ "4h.p.i",
                               str_detect(samples, "_12hrs") ~ "12h.p.i",
                               str_detect(samples, "_24hrs") ~ "24h.p.i",
                               str_detect(samples, "_72hrs") ~ "72h.p.i"),
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
        legend.position = c(0.9, 0.85),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 14, face = "bold"))

# plot transcript maps
plotTranscripts(unique(geneIDs(raw_data))[1], raw_data, 
                samples = c(pData(raw_data) %>% select(sample_name) %>% pull()),
                legend = T,
                labelTranscripts = T)


###### Gene counts with GenomicAlignments package
library(GenomicAlignments)
library(rtracklayer)

# read bam file
bam12s <- readGAlignments("results/hisat2/bulk/sortedTHEV_12hrsSamples.bam")

gtf_file <- rtracklayer::import("raw_files/annotations/thev_predicted_genes.gtf")

tallies <- summarizeOverlaps(gtf_file, bam12s1, mode = "Union",
                             ignore.strand = F, inter.feature = F,
                             singleEnd = F, fragments = T)






 
