#!/usr/bin/env zsh


filedir=raw_files/annotations

# remove older gff version

if [ -f $filedir/thev_predicted_genes.gff ] 
    then
        rm $filedir/thev_predicted_genes.gff
        echo "old GFF file removed"
    else
        echo "no GFF file in directory"
fi

################# CREATE .GFF3 FILE FROM .BED ###########
echo "Converting .bed file to .gff3 file ..." &&
agat_convert_bed2gff.pl --bed $filedir/THEVannotated_genesOnly.bed --source AbQuaye -o $filedir/thev_predicted_genes.gff && 
echo ".GFF3 file created successfully" || exit 1 