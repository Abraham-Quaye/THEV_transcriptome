#!/usr/bin/env zsh
set -e


zsh=scripts/zsh
r=scripts/r

#################### ACTIVATE CONDA ENVIRONMENT (RNASEQ) ##################
eval "$(conda shell.zsh hook)"
conda activate rnaseq && echo "successfull activation" || exit 1 


################# CREATE .GFF3 FILE FROM .BED ###########
$zsh/make_gff.zsh

################# CONVERT .GFF3 TO GTF ###########
$zsh/make_gtf.zsh

################### REMOVE UNWANTED GFF FEATURES ####################
R -e "source('$r/makingGFFfile.R')"

################### EXTRACT SPLICE-SITES ####################
$zsh/extract_ss.zsh

################### EXTRACT EXONS ####################
$zsh/extract_exons.zsh

#################### BUILD THEV GENOMIC INDEX FOR MAPPING WITH HISAT2 ######
$zsh/build_genome_index.zsh

#################### MOVE AGAT LOGFILE #################
rm raw_files/annotations/thev_predicted_genes.agat.log 
mv thev_predicted_genes.agat.log raw_files/annotations

#################### MAP READS TO THEV GENOME WITH HISAT2 #############
$zsh/map_sort_to_bam.zsh

#################### INDEX ALL SORTED .BAM FILES ##################
$zsh/index.zsh

#################### CONSTRUCT TRANSCRIPTS WITH STRINGTIE ##############
$zsh/assemble_transcripts.zsh 

#################### BULK MAP READS ##########
$zsh/bulk_map_sort_to_bam.zsh

#################### BULK COVERAGE #########
$zsh/bulk_coverage.zsh

#################### BULK DEPTH ############
$zsh/bulk_depth.zsh
# also removes the .bam files, leaving no bulk mapped files 

#################### MAKE FIGURES FOR DEPTH/COVERAGE ############
R -e "source('$r/thev_cov_depth.R')"

#################### DEACTIVATE ENVIRONMENT ##############
conda deactivate && echo "rna-seq env deactivated!!" 