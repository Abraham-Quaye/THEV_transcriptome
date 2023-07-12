#!/usr/bin/env zsh

mapdir=results/hisat2/bulk
covdir=results/hisat2/coverage

# Depth for 4hrs
(echo "Depth 4hrs starting..."
samtools depth -aH -r AY849321.1 $mapdir/sortedTHEV_4hrsSamples.bam > $covdir/thev_4hrsdepth.txt) & 


# Depth for 12hrs
(echo "Depth 12hrs starting..."
samtools depth -aH -r AY849321.1 $mapdir/sortedTHEV_12hrsSamples.bam > $covdir/thev_12hrsdepth.txt) &

# Depth for 24hrs
(echo "Depth 24hrs starting..."
samtools depth -aH -r AY849321.1 $mapdir/sortedTHEV_24hrsSamples.bam > $covdir/thev_24hrsdepth.txt) &

# Depth for 72hrs
(echo "Depth 72hrs starting..."
samtools depth -aH -r AY849321.1 $mapdir/sortedTHEV_72hrsSamples.bam > $covdir/thev_72hrsdepth.txt) &&

echo "All Depth estimations complete!!!" 
