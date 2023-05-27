#!/usr/bin/env zsh

mapdir=results/hisat2
fasta=raw_files/genome_file/AY849321.1.fa

echo "splice junction stats for 4h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_4S1.bed $mapdir/thev_sorted_4hrsS1.bam $fasta &
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_4S2.bed $mapdir/thev_sorted_4hrsS2.bam $fasta &
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_4S3.bed $mapdir/thev_sorted_4hrsS3.bam $fasta &

echo "splice junction stats for 12h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_12S1.bed $mapdir/thev_sorted_12hrsS1.bam $fasta &
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_12S3.bed $mapdir/thev_sorted_12hrsS3.bam $fasta &&

echo "splice junction stats for 24h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_24S1.bed $mapdir/thev_sorted_24hrsS1.bam $fasta &
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_24S2.bed $mapdir/thev_sorted_24hrsS2.bam $fasta &
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_24S3.bed $mapdir/thev_sorted_24hrsS3.bam $fasta &

echo "splice junction stats for 72h.p.i samples"
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_72S1.bed $mapdir/thev_sorted_72hrsS1.bam $fasta &
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_72S2.bed $mapdir/thev_sorted_72hrsS2.bam $fasta
regtools junctions extract -m 50 -s RF -o $mapdir/junction_stats_72S3.bed $mapdir/thev_sorted_72hrsS3.bam $fasta 

echo "splice junction stats complete!!!"
