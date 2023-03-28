
library(tidyverse)


# read primer data from all trancripts

primer_data <- list.files("wet_lab_validation/transcriptome_validation",
                          pattern = "[-a-zA-Z]+\\.txt",
                          full.names = TRUE)

all_primers <- map_dfr(primer_data, read_tsv, col_names = F, col_types = "c" )
colnames(all_primers) <- c(paste0("primer_", c("name", "seq", "length", "tm")))

## give the primer_tm column a standard format and filter out duplicates
all_primers <- all_primers %>%
  mutate(primer_tm = str_replace(string = primer_tm,
                                 pattern = "(\\w+\\s=\\s\\d+\\s\\w+)\\s\\(\\w+\\s?[\\w]+?\\)",
                                 replacement = "\\1")) %>% 
  unique()
  
## primers with unique names but identical sequences:
# all_primers[duplicated(all_primers$primer_seq),]
# 1. sMLP_trxptH_J2 R == truncMLP-trxptE R -> remove (sMLP_trxptH_J2 R) from primer list 

all_primers <- all_primers %>%
  filter(primer_name != "sMLP_trxptH_J2 R")

# Save Primer list to order:
export_table(all_primers, format = "html", sep = )















