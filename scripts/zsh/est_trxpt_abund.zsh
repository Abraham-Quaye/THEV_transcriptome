#!/usr/bin/env bash


mapdir=results/hisat2
GTF=results/gffcompare/updated_alltimes.combined.gtf
assembled=results/ballgown
samples="72hrsS1 72hrsS2 72hrsS3 24hrsS1 24hrsS2 24hrsS3 12hrsS1 12hrsS3 4hrsS1 4hrsS2"



echo "Estimating Abundances..."
for sample in $samples; do
    stringtie -p 10 -e -G $GTF -o $assembled/abund_${sample}/abund_${sample}.gtf $mapdir/thev_subset_${sample}.bam
done

echo "All Abundance Estimations Complete!"

