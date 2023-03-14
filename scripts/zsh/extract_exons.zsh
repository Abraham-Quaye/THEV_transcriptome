#!/bin/usr/env zsh


filedir=raw_files/annotations


if [ -f $filedir/thev_predicted_genes.exons ] 
    then
        rm $filedir/thev_predicted_genes.exons
        echo "old exon file removed"
    else
        echo "no exon file in directory"
fi


################# EXTRACT EXONS FROM GTF FILE 
echo "extracting exons ..."
extract_exons.py $filedir/thev_predicted_genes.gtf > $filedir/thev_predicted_genes.exons &&
echo "exon extration complete" || exit 4