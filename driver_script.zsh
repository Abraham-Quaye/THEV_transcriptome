#!/bin/zsh
set -e


scripts=scripts/zsh
#################### ACTIVATE CONDA ENVIRONMENT (RNASEQ) ##################
eval "$(conda shell.zsh hook)"
conda activate rnaseq && echo "successfull activation" || exit 1

#################### BUILD THEV GENOMIC INDEX FOR MAPPING WITH HISAT2 ######
$scripts/build_genome_index.zsh

: '
#################### MAP READS TO THEV GENOME WITH HISAT2 #############
$scripts/mapping.zsh

#################### SORT MAPPED READS AND CONVERT TO .BAM WITH SAMTOOLS ######
$scripts/sort_sam_to_bam.zsh

#################### INDEX ALL SORTED .BAM FILES ##################
$scripts/index.zsh

#################### DELETE ALL .SAM FILES ##################
$scripts/trash_sam.zsh

#################### CONSTRUCT TRANSCRIPTS WITH STRINGTIE ##############
$scripts/assemble_transcripts.zsh 

#################### BULK MAP READS ##########
$scripts/bulk_mapping.zsh

#################### BULK SORT .SAM TO .BAM ###########
$scripts/bulk_sort.zsh
# also deletes all .sam files, leaving only .bam files

#################### BULK COVERAGE #########
$scripts/bulk_coverage.zsh

#################### BULK DEPTH ############
$scripts/bulk_depth.zsh
# also removes the .bam files, leaving no bulk mapped files

conda deactivate && echo "rna-seq env deactivated!!" '