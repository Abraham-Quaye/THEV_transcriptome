#!/usr/bin/env zsh


mapdir=results/hisat2/bulk/uninfected

# Index for 4hrs
echo "Indexing Uninfected..."
samtools index $mapdir/sortedTHEV_72hrsNeg.bam &
samtools index $mapdir/sortedTHEV_24hrsNeg.bam ;
samtools index $mapdir/sortedTHEV_12hrsNeg.bam &
samtools index $mapdir/sortedTHEV_4hrsNeg.bam ;
echo "All Indexing complete!"
