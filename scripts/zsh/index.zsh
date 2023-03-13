#!/usr/bin/env zsh


mapdir=results/hisat2

########################################## SECTION III ##############################################
# indexing all bam files

# Index for 4hrs
echo "Indexing 4hrs..."
samtools index $mapdir/thev_sorted_4hrsS1.bam &
samtools index $mapdir/thev_sorted_4hrsS2.bam &
samtools index $mapdir/thev_sorted_4hrsS3.bam ;

# Index for 12hrs
echo "Indexing 12hrs..."
samtools index $mapdir/thev_sorted_12hrsS1.bam &
samtools index $mapdir/thev_sorted_12hrsS3.bam ;

# Index for 24hrs
echo "Index 24hrs starting..."
samtools index $mapdir/thev_sorted_24hrsS1.bam &
samtools index $mapdir/thev_sorted_24hrsS2.bam &
samtools index $mapdir/thev_sorted_24hrsS3.bam ;

# Index for 72hrs
echo "Index 72hrs starting..."
samtools index $mapdir/thev_sorted_72hrsS1.bam &
samtools index $mapdir/thev_sorted_72hrsS2.bam &
samtools index $mapdir/thev_sorted_72hrsS3.bam &&

echo "All Indexing complete!"
