#!/usr/bin/env zsh


mapdir=results/hisat2/bulk

# indexing all bam files

# Index for 4hrs
echo "Indexing 4hrs..."
samtools index $mapdir/sortedTHEV_4hrsSamples.bam &

# Index for 12hrs
echo "Indexing 12hrs..." &
samtools index $mapdir/sortedTHEV_12hrsSamples.bam &

# Index for 24hrs
echo "Index 24hrs starting..." &
samtools index $mapdir/sortedTHEV_24hrsSamples.bam &

# Index for 72hrs
echo "Index 72hrs starting..."&
samtools index $mapdir/sortedTHEV_72hrsSamples.bam ;

echo "All Indexing complete!"
