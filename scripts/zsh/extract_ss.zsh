#!/usr/bin/env zsh


filedir=raw_files/annotations


if [ -f $filedir/thev_predicted_genes.ss ] 
    then
        rm $filedir/thev_predicted_genes.ss
        echo "old splice-site file removed"
    else
        echo "no splice-site file in directory"
fi

################# EXTRACT SPLICESITES FROM GTF FILE ################
echo "extracting splice-sites ..."
extract_splice_sites.py $filedir/thev_predicted_genes.gtf > $filedir/thev_predicted_genes.ss &&
echo "splice-site extraction completed successfully" || exit 3