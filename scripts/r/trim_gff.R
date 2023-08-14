#!/usr/bin/env Rscript


library(tidyverse)


gff <-  read_tsv("raw_files/annotations/thev_predicted_genes.gff",
                 col_names = FALSE,
                 comment = "#",
                 show_col_types = FALSE) %>%
  filter(!X3 %in% c("three_prime_UTR", "five_prime_UTR")) %>% 
  separate_wider_delim(X9, delim = ";", names = c("id", "name", "others"),
                       too_many = "merge", too_few = "align_start") %>%
  separate_wider_delim(name, delim = "=", names = c("attr", "gene_name")) %>% 
  rename(chr = X1, source = X2, feature = X3) %>% 
  mutate(id = ifelse(feature == "gene", paste0("ID=", gene_name), id),
         gene_name = ifelse(attr == "Parent",
                            ifelse(lag(gene_name, n = 1) %in% c(as.character(c(1:23))),
                                   lag(gene_name, n = 2), lag(gene_name, n = 1)), 
                            gene_name),
         gene_name = ifelse(gene_name %in% c(as.character(c(1:23))), lag(gene_name, n = 1), gene_name),
         gene_name = ifelse(gene_name %in% c(as.character(c(1:23))), lag(gene_name, n = 1), gene_name),
         gene_attr = paste0("gene=", gene_name)
         ) %>%
  unite(col = "attr1", attr:gene_name, sep = "=", remove = T) %>% 
  unite(col = "attributes", id:gene_attr, sep = ";", na.rm = T)

# Save table
write.table(gff, "raw_files/annotations/thev_predicted_genes.gff",
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)





