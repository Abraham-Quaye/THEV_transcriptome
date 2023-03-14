#!/usr/bin/env zsh
set -e


zsh=scripts/zsh
r=scripts/r

#################### ACTIVATE CONDA ENVIRONMENT (RNASEQ) ##################
eval "$(conda shell.zsh hook)"
conda activate rnaseq && echo "successfull activation" || exit 1 


################# CREATE .GFF3 FILE FROM .BED ###########
$zsh/make_gff.zsh

################### REMOVE UNWANTED GFF FEATURES ####################
# R -e "source('$r/makingGFFfile.R')"

################# CONVERT .GFF3 TO GTF ###########
$zsh/make_gtf.zsh


#################### BUILD THEV GENOMIC INDEX FOR MAPPING WITH HISAT2 ######
$zsh/build_genome_index.zsh


#################### MAP READS TO THEV GENOME WITH HISAT2 #############
$zsh/mapping.zsh

#################### SORT MAPPED READS AND CONVERT TO .BAM WITH SAMTOOLS ######
$zsh/sort_sam_to_bam.zsh

#################### INDEX ALL SORTED .BAM FILES ##################
$zsh/index.zsh

#################### DELETE ALL .SAM FILES ##################
$zsh/trash_sam.zsh

#################### CONSTRUCT TRANSCRIPTS WITH STRINGTIE ##############
$zsh/assemble_transcripts.zsh 

#################### BULK MAP READS ##########
$zsh/bulk_mapping.zsh

#################### BULK SORT .SAM TO .BAM ###########
$zsh/bulk_sort.zsh
# also deletes all .sam files, leaving only .bam files

#################### BULK COVERAGE #########
$zsh/bulk_coverage.zsh

#################### BULK DEPTH ############
$zsh/bulk_depth.zsh
# also removes the .bam files, leaving no bulk mapped files 

#################### MAKE FIGURES FOR DEPTH/COVERAGE ############
R -e "source('$r/thev_cov_depth.R')"

#################### DEACTIVATE ENVIRONMENT ##############
conda deactivate && echo "rna-seq env deactivated!!" 