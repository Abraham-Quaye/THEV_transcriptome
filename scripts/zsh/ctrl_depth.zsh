#!/usr/bin/env zsh

mapdir=results/hisat2/bulk/uninfected
covdir=results/hisat2/coverage

# Depth for 4hrs
(echo "Depth 4hrs starting..."
samtools depth -aH $mapdir/sortedTHEV_4hrsNeg.bam > $covdir/ctrl_4hrsdepth.txt) & 


# Depth for 12hrs
(echo "Depth 12hrs starting..."
samtools depth -aH $mapdir/sortedTHEV_12hrsNeg.bam > $covdir/ctrl_12hrsdepth.txt) &

# Depth for 24hrs
(echo "Depth 24hrs starting..."
samtools depth -aH $mapdir/sortedTHEV_24hrsNeg.bam > $covdir/ctrl_24hrsdepth.txt) &

# Depth for 72hrs
(echo "Depth 72hrs starting..."
samtools depth -aH $mapdir/sortedTHEV_72hrsNeg.bam > $covdir/ctrl_72hrsdepth.txt) &&

echo "All Depth estimations complete!!!" 
