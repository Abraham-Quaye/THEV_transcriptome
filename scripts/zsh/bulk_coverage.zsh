#!/usr/bin/env zsh


mapdir=results/hisat2/bulk
covfile=results/hisat2/coverage/bulk_coverage.txt

# coverage for 4hrs
echo "Coverage 4hrs starting..." &&
echo "Coverage 4hrs" > $covfile &
samtools coverage $mapdir/sortedTHEV_4hrsSamples.bam >> $covfile &&
echo "Coverage 4hrs complete!" &

# coverage for 12hrs
echo "Coverage 12hrs starting..." &
echo "Coverage 12hrs" >> $covfile &
samtools coverage $mapdir/sortedTHEV_12hrsSamples.bam >> $covfile &&
echo "Coverage 12hrs complete!" &

# coverage for 24hrs
echo "Coverage 24hrs starting..." &
echo "Coverage 24hrs" >> $covfile &
samtools coverage $mapdir/sortedTHEV_24hrsSamples.bam >> $covfile &&
echo "Coverage 24hrs complete!" &

# coverage for 72hrs
echo "Coverage 72hrs starting..." &
echo "Coverage 72hrs" >> $covfile &
samtools coverage $mapdir/sortedTHEV_72hrsSamples.bam >> $covfile &&
echo "All coverage estimation complete!!!"
