#!/usr/bin/env zsh


mapdir=results/hisat2/bulk

# SECTION II -- BULK-SORTING

# alert for start
#sort 4hrs sample
echo "Sorting 4hrs ..."
samtools sort -@ 10  -o $mapdir/sortedTHEV_4hrsSamples.bam $mapdir/thev_4hrsSamples.sam &

#sort 12hrs sample
echo "Sorting 12hrs..."
samtools sort -@ 10  -o $mapdir/sortedTHEV_12hrsSamples.bam $mapdir/thev_12hrsSamples.sam ;

#sort 24hrs sample
echo "Sorting 24hrs..."
samtools sort -@ 10  -o $mapdir/sortedTHEV_24hrsSamples.bam $mapdir/thev_24hrsSamples.sam &

#sort 72hrs sample
echo "Sorting 72hrs..."
samtools sort -@ 10  -o $mapdir/sortedTHEV_72hrsSamples.bam $mapdir/thev_72hrsSamples.sam ;

# alert for end
echo "Sorting Complete!"


# SECTION III -- TRASH .SAM FILES
rm $mapdir/*.sam
