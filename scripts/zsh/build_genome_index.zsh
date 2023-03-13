#!/bin/zsh


filedir=raw_files/annotations
genomedir=raw_files/genome_file
idxdir=raw_files/thevgenome_index


################# CREATE .GFF3 FILE FROM .BED ###########
echo "Converting .bed file to .gff3 file ..." &&
agat_convert_bed2gff.pl --bed $filedir/THEVannotated_genesOnly.bed --source abquaye --primary_tag CDS -o $filedir/thev_predicted_genes.gff && 
echo ".GFF3 file created successfully" || exit 1

################# CONVERT .GFF3 TO GTF ###########
echo "Converting .gff3 file to .gtf file ..." ;
agat_convert_sp_gff2gtf.pl --gff $filedir/thev_predicted_genes.gff -o $filedir/thev_predicted_genes.gtf && 
echo ".GTF file created successfully" || exit 2


################# EXTRACT SPLICESITES FROM GTF FILE ################
echo "extracting splice-sites ..."
extract_splice_sites.py $filedir/thev_predicted_genes.gff > $filedir/thev_predicted_genes.ss &&
echo "splice-site extraction completed successfully" || exit 3


################# EXTRACT EXONS FROM GTF FILE 
echo "extracting exons ..."
extract_exons.py $filedir/thev_predicted_genes.gff > $filedir/thev_predicted_genes.exons &&
echo "exon extration complete" || exit 4


################# COMMAND TO BUILD THEV GENOMIC INDEX #############
echo "Building thev genomic index ..."
hisat2-build -p 10 --ss $filedir/thev_predicted_genes.ss --exon $filedir/thev_predicted_genes.exons $genomedir/THEV.fa $idxdir/thev_tran &&
echo "Index built successfully" || exit 5 
