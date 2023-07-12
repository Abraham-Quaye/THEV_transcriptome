#!/usr/bin/env zsh


mapdir=results/hisat2/bulk

# indexing all bam files

# Index for 4hrs
echo "Indexing 4hrs..."
samtools index $mapdir/subsetTHEV_4hrs.bam &

# Index for 12hrs
echo "Indexing 12hrs..." &
samtools index $mapdir/subsetTHEV_12hrs.bam &

# Index for 24hrs
echo "Index 24hrs starting..." &
samtools index $mapdir/subsetTHEV_24hrs.bam &

# Index for 72hrs
echo "Index 72hrs starting..."&
samtools index $mapdir/subsetTHEV_72hrs.bam ;

echo "All Indexing complete!"
