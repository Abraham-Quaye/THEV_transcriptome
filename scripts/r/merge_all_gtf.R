######## Merge gtf files of all time-points from stringtie ------

library(tidyverse)

all_gtfs <- list.files(path = "results/stringtie",
                       pattern = "[a-z_]+\\d+[a-zA-Z]+\\d\\.gtf",
                       full.names = TRUE)

names(all_gtfs) <- c("t12s1", "t12s2", "t24s1", "t24s2", "t24s3",
                     "t4s1", "t4s2", "t4s3", "t72s1", "t72s2", "t72s3")

all_merged <- map_dfr(all_gtfs, read_tsv, comment = "#", col_names = FALSE)

write.table(all_merged, "results/stringtie/all_real_transcripts_mergedRR.gtf",
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)



