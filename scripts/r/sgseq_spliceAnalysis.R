
library(SGSeq)
library(tidyverse)

## make list of samples with corresponding bam files to extract information
# sample <- tibble(sample_name = c(paste0("12hrs", c("S1", "S3")),
#                                  paste0(c("24hrs"), c("S1", "S2", "S3")),
#                                  paste0(c("4hrs"), c("S1", "S2", "S3")),
#                                  paste0(c("72hrs"), c("S1", "S2", "S3"))),
#                  file_bam = list.files("results/hisat2", "^thev_sorted_.+bam$",
#                                        full.names = T)
#                  )

# ## extract information from bam files. This is a slow function so I'll save the data
# # as a tsv and load it henceforth
# sample_info <- getBamInfo(sample, cores = 8)
# write.table(sample_info, "results/sgseq/bam_extracted_sampleData.tsv",
#             col.names = T, row.names = F, quote = F, sep = "\t")

### ====================================
# load in sample/bam information

bam_info <- read_tsv("results/sgseq/bam_extracted_sampleData.tsv", show_col_types = F) %>%
  slice(3)

### ====================================
# load in predicted genes from ncbi

features <- importTranscripts("raw_files/annotations/thev_from_NCBI.gtf")
dtx_features <- convertToTxFeatures(features)
sg_features <- convertToSGFeatures(tx_features)

### ====================================
# analyze splice sites

spliceSites24 <- analyzeFeatures(bam_info, tx_features, cores = 8, verbose = T)


# visualize
plotFeatures(spliceSites24, geneID = 1, heightPanels = c(3, 1), square = T,
             toscale = "exon",
             col = colorRampPalette(c("blue", "red"))(256), cex = 0.5)




