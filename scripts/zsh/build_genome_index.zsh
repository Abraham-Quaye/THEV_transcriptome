#!/usr/bin/env zsh


filedir=raw_files/annotations
genomedir=raw_files/genome_file
idxdir=raw_files/thevgenome_index

# remove older files

if [ -f $idxdir/thev_tran.7.ht2 ] 
    then
        rm $idxdir/*.ht2
        echo "old index files removed"
    else
        echo "no index file in directory"
fi


################# COMMAND TO BUILD THEV GENOMIC INDEX #############
echo "Building thev genomic index ..."
hisat2-build -p 10 --ss $filedir/thev_predicted_genes.ss --exon $filedir/thev_predicted_genes.exons $genomedir/AY849321.1.fa $idxdir/thev_tran &&
echo "Index built successfully" || exit 5 