#!/usr/bin/env zsh


mapdir=results/hisat2/bulk

# indexing all bam files

# Index for 4hrs
echo "Subsetting 4hrs..."
samtools view -b -@ 8 $mapdir/sortedTHEV_4hrsSamples.bam AY849321.1 > $mapdir/subsetTHEV_4hrs.bam

# Index for 12hrs
echo "Subsetting 12hrs..."
samtools view -b -@ 8 $mapdir/sortedTHEV_12hrsSamples.bam AY849321.1 > $mapdir/subsetTHEV_12hrs.bam

# Index for 24hrs
echo "Subsetting 24hrs..."
samtools view -b -@ 8 $mapdir/sortedTHEV_24hrsSamples.bam AY849321.1 > $mapdir/subsetTHEV_24hrs.bam

# Index for 72hrs
echo "Subsetting 72hrs..."
samtools view -b -@ 8 $mapdir/sortedTHEV_72hrsSamples.bam AY849321.1 > $mapdir/subsetTHEV_72hrs.bam

echo "All Subsetting complete!"
