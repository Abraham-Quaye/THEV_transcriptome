#!/usr/bin/env zsh

mapdir=results/hisat2/bulk
fasta=raw_files/genome_file/AY849321.1.fa

echo "splice junction stats for 4h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_4hrs.bed $mapdir/subsetTHEV_4hrs.bam $fasta &

echo "splice junction stats for 12h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_12hrs.bed $mapdir/subsetTHEV_12hrs.bam $fasta ;

echo "splice junction stats for 24h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_24hrs.bed $mapdir/subsetTHEV_24hrs.bam $fasta &

echo "splice junction stats for 72h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_72hrs.bed $mapdir/subsetTHEV_72hrs.bam $fasta ;

echo "Bulk splice junction stats complete!!!"
