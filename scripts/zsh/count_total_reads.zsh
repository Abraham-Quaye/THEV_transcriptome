#!/usr/bin/env zsh

countfile=results/hisat2/coverage/bulk_counts.txt
mapdir=results/hisat2/bulk

echo "counting total reads"
echo "total_reads" > $countfile ;
samtools view -c $mapdir/sortedTHEV_4hrsSamples.bam >> $countfile ;

samtools view -c $mapdir/sortedTHEV_12hrsSamples.bam >> $countfile ;

samtools view -c $mapdir/sortedTHEV_24hrsSamples.bam >> $countfile ;

samtools view -c $mapdir/sortedTHEV_72hrsSamples.bam >> $countfile ;

echo "counting reads complete"