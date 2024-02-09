#!/usr/bin/env Rscript

library(magrittr)
library(scales)
library(plyr)
library(rtracklayer)
library(glue)
library(devtools)
library(genefilter)
library(tidyverse)


# ====================================================================
# TRANSCRIPT ABUNDANCE ANALYSIS
# ====================================================================
# load in abund files
abund_metadata <- tibble(dir = c(list.dirs("results/abundances", full.names = T) %>% .[-1]),
                     abund_file = map_chr(dir, function(dir) list.files(path = dir,
                                             pattern = "^abund_\\d{1,2}hrsS\\d\\.gtf$",
                                             full.names = T)),
                     raw_tabs = map(abund_file, \(abund_file) as_tibble(import(abund_file))),
                     time_labs = case_when(str_detect(dir, "_4hrs") ~ "4hpi",
                                           str_detect(dir, "_12hrs") ~ "12hpi",
                                           str_detect(dir, "_24hrs") ~ "24hpi",
                                           str_detect(dir, "_72hrs") ~ "72hpi",
                                           TRUE ~ NA_character_)
                     )

# join the dataframes as one
abund_data <- abund_metadata$raw_tabs %>%
  set_names(abund_metadata$time_labs) %>%
  bind_rows(., .id = "timepoint") %>%
  rename(fpkm = FPKM, tpm = TPM) %>%
  filter(type == "transcript") %>%
  mutate(tpm = as.numeric(tpm),
         fpkm = as.numeric(fpkm)) %>%
  group_by(timepoint, transcript_id, start, end) %>%
  reframe(total_tmp = sum(tpm)) %>%
  arrange(start)
  

abund_data %>%
  ggplot(aes(transcript_id, total_tmp, fill = timepoint)) +
  geom_col(position = "dodge") +
  theme(axis.text.x = element_text(angle = 45))