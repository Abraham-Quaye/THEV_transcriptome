#!/usr/bin/env zsh

mapdir=results/hisat2/bulk
fasta=raw_files/genome_file/AY849321.1.fa

echo "splice junction stats for 4h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_4hrs.bed $mapdir/sortedTHEV_4hrsSamples.bam $fasta &

echo "splice junction stats for 12h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_12hrs.bed $mapdir/sortedTHEV_12hrsSamples.bam $fasta ;

echo "splice junction stats for 24h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_24hrs.bed $mapdir/sortedTHEV_24hrsSamples.bam $fasta &

echo "splice junction stats for 72h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_72hrs.bed $mapdir/sortedTHEV_72hrsSamples.bam $fasta ;

echo "Bulk splice junction stats complete!!!"
