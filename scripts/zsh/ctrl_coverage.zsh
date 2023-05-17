#!/usr/bin/env zsh


mapdir=results/hisat2/bulk/uninfected
covfile=results/hisat2/coverage/ctrl_coverage.tsv

# coverage for 4hrs
echo "Coverage 4hrs starting..." &&
echo "Coverage 4hrs" > $covfile &
samtools coverage $mapdir/sortedTHEV_4hrsNeg.bam >> $covfile &&
echo "Coverage 4hrs complete!" &

# coverage for 12hrs
echo "Coverage 12hrs starting..." &
echo "Coverage 12hrs" >> $covfile &
samtools coverage $mapdir/sortedTHEV_12hrsNeg.bam >> $covfile &&
echo "Coverage 12hrs complete!" &

# coverage for 24hrs
echo "Coverage 24hrs starting..." &
echo "Coverage 24hrs" >> $covfile &
samtools coverage $mapdir/sortedTHEV_24hrsNeg.bam >> $covfile &&
echo "Coverage 24hrs complete!" &

# coverage for 72hrs
echo "Coverage 72hrs starting..." &
echo "Coverage 72hrs" >> $covfile &
samtools coverage $mapdir/sortedTHEV_72hrsNeg.bam >> $covfile &&
echo "All coverage estimation complete!!!"
