#!/usr/bin/env zsh


mapdir=results/hisat2

if [ -f $mapdir/thev_sorted_4hrsS2.bam ] 
    then
        rm $mapdir/*.bam
        rm $mapdir/*.bai
        echo "old .bam and .bai files removed"
    else
        echo "no .bam files in directory"
fi
#################################### SECTION II #########################################
# sorting all mapped reads and converting to .bam files with samtools

#sort 4hrs sample
echo "Sorting .sam files and converting to .bam ..."
samtools sort -@ 8 -o $mapdir/thev_sorted_4hrsS1.bam $mapdir/thev_4hrsS1.sam &
samtools sort -@ 8 -o $mapdir/thev_sorted_4hrsS2.bam $mapdir/thev_4hrsS2.sam &
samtools sort -@ 8 -o $mapdir/thev_sorted_4hrsS3.bam $mapdir/thev_4hrsS3.sam ; 
echo "4hrs sorting done" ;

#sort 12hrs sample
samtools sort -@ 8 -o $mapdir/thev_sorted_12hrsS1.bam $mapdir/thev_12hrsS1.sam &
samtools sort -@ 8 -o $mapdir/thev_sorted_12hrsS3.bam $mapdir/thev_12hrsS3.sam &&
echo "12hrs sorting done" ;

#sort 24hrs sample
samtools sort -@ 8 -o $mapdir/thev_sorted_24hrsS1.bam $mapdir/thev_24hrsS1.sam &
samtools sort -@ 8 -o $mapdir/thev_sorted_24hrsS2.bam $mapdir/thev_24hrsS2.sam &
samtools sort -@ 8 -o $mapdir/thev_sorted_24hrsS3.bam $mapdir/thev_24hrsS3.sam ;
echo "24hrs sorting complete" ;

#sort 72hrs sample
samtools sort -@ 8 -o $mapdir/thev_sorted_72hrsS1.bam $mapdir/thev_72hrsS1.sam &
samtools sort -@ 8 -o $mapdir/thev_sorted_72hrsS2.bam $mapdir/thev_72hrsS2.sam &
samtools sort -@ 8 -o $mapdir/thev_sorted_72hrsS3.bam $mapdir/thev_72hrsS3.sam &&
echo "All sorting Complete!"