#!/usr/bin/env zsh


filedir=raw_files/annotations


################# CREATE .GFF3 FILE FROM .BED ###########
echo "Converting .bed file to .gff3 file ..." &&
agat_convert_bed2gff.pl --bed $filedir/THEVannotated_genesOnly.bed --source Genbank -o $filedir/thev_predicted_genes.gff && 
echo ".GFF3 file created successfully"

