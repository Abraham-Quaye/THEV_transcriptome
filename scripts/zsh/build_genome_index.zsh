#!/usr/bin/env zsh


filedir=raw_files/annotations
genomedir=raw_files/genome_file
idxdir=raw_files/thevgenome_index

################# EXTRACT SPLICESITES FROM GTF FILE ################
echo "extracting splice-sites ..."
extract_splice_sites.py $filedir/thev_predicted_genes.gtf > $filedir/thev_predicted_genes.ss &&
echo "splice-site extraction completed successfully" || exit 3


################# EXTRACT EXONS FROM GTF FILE 
echo "extracting exons ..."
extract_exons.py $filedir/thev_predicted_genes.gtf > $filedir/thev_predicted_genes.exons &&
echo "exon extration complete" || exit 4


################# COMMAND TO BUILD THEV GENOMIC INDEX #############
echo "Building thev genomic index ..."
hisat2-build -p 10 --ss $filedir/thev_predicted_genes.ss --exon $filedir/thev_predicted_genes.exons $genomedir/THEV.fa $idxdir/thev_tran &&
echo "Index built successfully" || exit 5 
