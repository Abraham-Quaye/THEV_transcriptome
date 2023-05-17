#!/usr/bin/env zsh

gffcompare -N -G -r raw_files/annotations/thev_from_NCBI.gtf \
-o results/stringtie/gffcomp_alltimes \
{results/stringtie/thev_24hrsS1.gtf,results/stringtie/thev_24hrsS2.gtf,\
results/stringtie/thev_24hrsS3.gtf,results/stringtie/thev_4hrsS2.gtf,\
results/stringtie/thev_4hrsS1.gtf,results/stringtie/thev_4hrsS3.gtf,\
results/stringtie/thev_12hrsS1.gtf,results/stringtie/thev_12hrsS3.gtf,\
results/stringtie/thev_72hrsS1.gtf,results/stringtie/thev_72hrsS2.gtf,\
results/stringtie/thev_72hrsS3.gtf}