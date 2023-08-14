#!/usr/bin/env Rscript


library(tidyverse)


gtf <-  read_tsv("raw_files/annotations/thev_predicted_genes.gtf",
                 col_names = FALSE,
                 comment = "#",
                 show_col_types = FALSE) %>%
  filter(X3 == "gene") %>% 
  select(X1, X4, X5, X7) %>% 
  mutate(X4 = X4 -1,
         X5 = X5 -1)


write.table(gtf, "raw_files/annotations/thev_predicted_genes.exons",
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE
)





