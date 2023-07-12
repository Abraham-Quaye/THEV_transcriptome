#!/usr/bin/env zsh


mapdir=results/hisat2

# Subset for 4hrs
echo "Subsetting 4hrs..."
samtools view -b -@ 8 $mapdir/thev_sorted_4hrsS1.bam AY849321.1 > $mapdir/thev_subset_4hrsS1.bam
samtools view -b -@ 8 $mapdir/thev_sorted_4hrsS2.bam AY849321.1 > $mapdir/thev_subset_4hrsS2.bam
samtools view -b -@ 8 $mapdir/thev_sorted_4hrsS3.bam AY849321.1 > $mapdir/thev_subset_4hrsS3.bam

# Subset for 12hrs
echo "Subsetting 12hrs..."
samtools view -b -@ 8 $mapdir/thev_sorted_12hrsS1.bam AY849321.1 > $mapdir/thev_subset_12hrsS1.bam
samtools view -b -@ 8 $mapdir/thev_sorted_12hrsS3.bam AY849321.1 > $mapdir/thev_subset_12hrsS3.bam

# Subset for 24hrs
echo "Subset 24hrs starting..."
samtools view -b -@ 8 $mapdir/thev_sorted_24hrsS1.bam AY849321.1 > $mapdir/thev_subset_24hrsS1.bam
samtools view -b -@ 8 $mapdir/thev_sorted_24hrsS2.bam AY849321.1 > $mapdir/thev_subset_24hrsS2.bam
samtools view -b -@ 8 $mapdir/thev_sorted_24hrsS3.bam AY849321.1 > $mapdir/thev_subset_24hrsS3.bam

# Subset for 72hrs
echo "Subset 72hrs starting..."&
samtools view -b -@ 8 $mapdir/thev_sorted_72hrsS1.bam AY849321.1 > $mapdir/thev_subset_72hrsS1.bam
samtools view -b -@ 8 $mapdir/thev_sorted_72hrsS2.bam AY849321.1 > $mapdir/thev_subset_72hrsS2.bam
samtools view -b -@ 8 $mapdir/thev_sorted_72hrsS3.bam AY849321.1 > $mapdir/thev_subset_72hrsS3.bam

echo "All Subsetting complete!"
