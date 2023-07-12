#!/usr/bin/env zsh


filedir=raw_files/annotations

################# CONVERT THEV .GFF3 TO GTF ###########
echo "Converting THEV .gff3 file to .gtf file ..." ;
agat_convert_sp_gff2gtf.pl --gff $filedir/thev_from_NCBI.gff3 -o $filedir/thev_from_NCBI.gtf && 
echo "THEV .GTF file created successfully"


################# CONVERT TURKEY .GFF3 TO GTF ###########
echo "Converting TURKEY .gff file to .gtf file ..." ;
agat_convert_sp_gff2gtf.pl --gff $filedir/turkey_genome.gff -o $filedir/turkey_genome.gtf ;
echo "TURKEY .GTF file created successfully" ;


mv thev_from_NCBI.agat.log $filedir ;
mv turkey_genome.agat.log $filedir ;
