#!/usr/bin/env zsh
dir=results/stringtie
gffcompare -G -r raw_files/annotations/thev_from_NCBI.gtf \
-o results/gffcompare/gffcomp_alltimes \
{$dir/transcripts_merged_4hrs.gtf,$dir/transcripts_merged_12hrs.gtf,\
$dir/transcripts_merged_24hrs.gtf,$dir/transcripts_merged_72hrs.gtf}

echo "GFFCOMPARE complete"