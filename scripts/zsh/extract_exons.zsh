#!/usr/bin/env zsh


filedir=raw_files/annotations



################# EXTRACT EXONS FROM THEV GTF FILE 
echo "THEV extracting exons ..."
extract_exons.py $filedir/thev_from_NCBI.gtf > $filedir/thev_predicted_genes.exons &&
echo "THEV exon extration complete"

################# EXTRACT EXONS FROM TURKEY GTF FILE 
echo "TURKEY extracting exons ..."
extract_exons.py $filedir/turkey_genome.gtf > $filedir/turkey_genome.exons &&
echo "TURKEY exon extration complete"

cat $filedir/*.exons > $filedir/host_thev.exons