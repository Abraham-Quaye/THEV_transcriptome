#!/usr/bin/env zsh


mapdir=results/hisat2/bulk
covdir=results/hisat2/coverage

# coverage for 4hrs
(echo "Coverage 4hrs starting..." &&
echo "Coverage 4hrs" > $covdir/host_thev_cov4.txt &
samtools coverage $mapdir/sortedTHEV_4hrsSamples.bam >> $covdir/host_thev_cov4.txt) &

# coverage for 12hrs
(echo "Coverage 12hrs starting..." &
echo "Coverage 12hrs" > $covdir/host_thev_cov12.txt &
samtools coverage $mapdir/sortedTHEV_12hrsSamples.bam >> $covdir/host_thev_cov12.txt) &

# coverage for 24hrs
(echo "Coverage 24hrs starting..." &
echo "Coverage 24hrs" > $covdir/host_thev_cov24.txt
samtools coverage $mapdir/sortedTHEV_24hrsSamples.bam >> $covdir/host_thev_cov24.txt) &

# coverage for 72hrs
(echo "Coverage 72hrs starting..." &
echo "Coverage 72hrs" > $covdir/host_thev_cov72.txt &
samtools coverage $mapdir/sortedTHEV_72hrsSamples.bam >> $covdir/host_thev_cov72.txt) ;

echo "All coverage estimation complete!!!"
