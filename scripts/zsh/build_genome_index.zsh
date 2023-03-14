#!/usr/bin/env zsh


filedir=raw_files/annotations
genomedir=raw_files/genome_file
idxdir=raw_files/thevgenome_index

# remove older files

if [ -f $filedir/thev_predicted_genes.ss ] 
    then
        rm $filedir/thev_predicted_genes.ss
        echo "old splice-site file removed"
    else
        echo "no splice-site file in directory"
fi
if [ -f $filedir/thev_predicted_genes.exons ] 
    then
        rm $filedir/thev_predicted_genes.exons
        echo "old exon file removed"
    else
        echo "no exon file in directory"
fi
if [ -f $idxdir/thev_tran.7.ht2 ] 
    then
        rm $idxdir/*.ht2
        echo "old index files removed"
    else
        echo "no index file in directory"
fi


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