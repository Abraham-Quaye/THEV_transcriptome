library(tidyverse)

# gene <- read_tsv("THEVgenelist_uniprot.tsv")
# 
# anno <- read_tsv("THEVannotated_genesOnly.txt", col_names = F, col_types = cols(X11 = "c"))
# 
# #write.table(anno, "THEVannotated_genesOnly.bed", row.names = F, col.names = F, quote = F)
# 
# anno <- rename(anno, chr = X1, chrStart = X2, chrEnd = X3, gene_name = X4, score = X5, strand = X6, thickStart = X7, thickEnd = X8,
#                itemRGB = X9, blockCounts = X10, blockSizes = X11, blockStarts = X12)
# 
# 
# basket <- c()
# for(i in 1:nrow(anno)){
#   basket[i] <- (anno$thickStart[i] + 1)
# }
# 
# basket2 <- c()
# for(i in 1:nrow(anno)){
#   basket2[i] <- (anno$thickEnd[i] + 1)
# }
# 
# newgff <- data.frame(seqname = anno$gene_name, source = rep("UniprotKB", nrow(anno)), feature = rep("CDS", nrow(anno)),
#                      start = basket, end = basket2, score = rep(".", nrow(anno)), strand = anno$strand, phase = rep(0,nrow(anno)))
# 
# #write.table(newgff, "thevgenes.txt", quote = F, row.names = F, sep = "\t")
# 
# pregff <- read_tsv("thevgenes.txt")
# write.table(pregff, "thevgenes.gff", quote = F, row.names = F, col.names = F,sep = "\t")
# 
# 
# # exons <- read_tsv("Book2.txt")
# # 
# # write.table(exons, "thevexons.gtf", quote = F, row.names = F, sep = "\t")
# 
# 
bed <- read_tsv("THEVannotated_genesOnly.txt",  col_names = F, col_types = cols(X11 = "c"))
write.table(bed, "THEVannotated_genesOnly.bed", col.names = F, row.names = F, sep = "\t", quote = F)
# 
# bed2 <- read_tsv("THEVannotated_genesOnly.bed", col_names = F, col_types = cols(X1d1 = "c"))





