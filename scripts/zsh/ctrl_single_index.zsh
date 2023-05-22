#!/usr/bin/env zsh


mapdir=results/hisat2/map_mock

########################################## SECTION III ##############################################
# indexing all bam files

# Index for 4hrs
echo "Indexing ctrl 4hrs..."
samtools index $mapdir/ctrl_sorted_4N1.bam &
samtools index $mapdir/ctrl_sorted_4N2.bam &

# Index for 12hrs
echo "Indexing ctrl 12hrs..." &
samtools index $mapdir/ctrl_sorted_12N1.bam &
samtools index $mapdir/ctrl_sorted_12N2.bam &

# Index for 24hrs
echo "Index ctrl 24hrs starting..." &
samtools index $mapdir/ctrl_sorted_24N1.bam &
samtools index $mapdir/ctrl_sorted_24N2.bam &

# Index for 72hrs
echo "Index ctrl 72hrs starting..."&
samtools index $mapdir/ctrl_sorted_72N1.bam &
samtools index $mapdir/ctrl_sorted_72N2.bam &&

echo "All ctrl indexing complete!"
