#!/usr/bin/env zsh


filedir=raw_files/annotations

# remove older gtf version

if [ -f $filedir/thev_from_NCBI.gtf ] 
    then
        rm $filedir/thev_predicted_genes.gtf
        echo "old GTF file removed"
    else
        echo "no GTF file in directory"
fi

################# CONVERT .GFF3 TO GTF ###########
echo "Converting .gff3 file to .gtf file ..." ;
agat_convert_sp_gff2gtf.pl --gff $filedir/thev_from_NCBI.gff3 -o $filedir/thev_from_NCBI.gtf && 
echo ".GTF file created successfully" || exit 2 


if [ -f $filedir/thev_predicted_genes.agat.log ]
    then
        rm $filedir/thev_from_NCBI.agat.log 
        echo "Old Agat logfile removed"
    else:
        echo "Agat logfile not yet moved"
fi

mv thev_from_NCBI.agat.log $filedir