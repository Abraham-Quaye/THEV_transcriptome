#!/bin/zsh

mapdir=results/hisat2/bulk
covdir=results/hisat2/coverage

# Depth for 4hrs
echo "Depth 4hrs starting..."
samtools depth -aH $mapdir/sortedNCBI_4hrsSamples.bam > $covdir/mappedNCBI_4hrsdepth.txt
echo "Depth 4hrs complete!"

# Depth for 12hrs
echo "Depth 12hrs starting..."
samtools depth -aH $mapdir/sortedNCBI_12hrsSamples.bam > $covdir/mappedNCBI_12hrsdepth.txt
echo "Depth 12hrs complete!"

# Depth for 24hrs
echo "Depth 24hrs starting..."
samtools depth -aH $mapdir/sortedNCBI_24hrsSamples.bam > $covdir/mappedNCBI_24hrsdepth.txt
echo "Depth 24hrs complete!"

# Depth for 72hrs
echo "Depth 72hrs starting..."
samtools depth -aH $mapdir/sortedNCBI_72hrsSamples.bam > $covdir/mappedNCBI_72hrsdepth.txt
echo "All Depth estimations complete!!!"


# SECTION III -- TRASH BULK .BAM FILES
rm $mapdir/*.bam