#!/usr/bin/env zsh


filedir=raw_files/annotations
genomedir=raw_files/genome_file
idxdir=raw_files/host_virus_genome_index


################# COMMAND TO BUILD THEV GENOMIC INDEX #############
# echo "Building thev genomic index ..."
# hisat2-build -p 10 --ss $filedir/thev_predicted_genes.ss --exon $filedir/thev_predicted_genes.exons $genomedir/AY849321.1.fa $idxdir/thev_tran &&
# echo "Index built successfully" || exit 5 


################# COMMAND TO BUILD THEV GENOMIC INDEX #############
echo "Building THEV AND TURKEY genomic index ..."
hisat2-build -p 10 --ss $filedir/host_thev.ss --exon $filedir/host_thev.exons $genomedir/mix_genomes.fa $idxdir/mix_tran ;


# $genomedir/mp_chr1.fa,$genomedir/mp_chr2.fa,$genomedir/mp_chr3.fa,$genomedir/mp_chr4.fa,$genomedir/mp_chr5.fa,$genomedir/mp_chr6.fa,$genomedir/mp_chr7.fa,$genomedir/mp_chr8.fa,$genomedir/mp_chr9.fa,$genomedir/mp_chr10.fa,$genomedir/mp_chr11.fa,$genomedir/mp_chr12.fa,$genomedir/mp_chr13.fa,$genomedir/mp_chr14.fa,$genomedir/mp_chr15.fa,$genomedir/mp_chr16.fa,$genomedir/mp_chr17.fa,$genomedir/mp_chr18.fa,$genomedir/mp_chr19.fa,$genomedir/mp_chr20.fa,$genomedir/mp_chr21.fa,$genomedir/mp_chr22.fa,$genomedir/mp_chr23.fa,$genomedir/mp_chr24.fa,$genomedir/mp_chr25.fa,$genomedir/mp_chr26.fa,$genomedir/mp_chr27.fa,$genomedir/mp_chr28.fa,$genomedir/mp_chr29.fa,$genomedir/mp_chr30.fa,$genomedir/mp_chrMT.fa,$genomedir/mp_chrW.fa,$genomedir/mp_chrZ.fa, 

echo "Index built successfully"
