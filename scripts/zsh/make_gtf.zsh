#!/usr/bin/env zsh


filedir=raw_files/annotations

# remove older gtf version

if [ -f $filedir/thev_predicted_genes.gtf ] 
    then
        rm $filedir/thev_predicted_genes.gtf
        echo "old GTF file removed"
    else
        echo "no GTF file in directory"
fi

################# CONVERT .GFF3 TO GTF ###########
echo "Converting .gff3 file to .gtf file ..." ;
agat_convert_sp_gff2gtf.pl --gff $filedir/thev_predicted_genes.gff -o $filedir/thev_predicted_genes.gtf && 
echo ".GTF file created successfully" || exit 2 