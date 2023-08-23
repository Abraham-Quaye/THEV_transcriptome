#!/usr/bin/env zsh

fasta=raw_files/genome_file/AY849321.1.fa
outdir=results/hisat2/bulk
gtf=results/gffcompare/updated_alltimes.combined.gtf

regtools junctions annotate -S -o $outdir/annot_4hrsSS.txt $outdir/junction_stats_4hrs.bed $fasta $gtf &

regtools junctions annotate -S -o $outdir/annot_12hrsSS.txt $outdir/junction_stats_12hrs.bed $fasta $gtf ;

regtools junctions annotate -S -o $outdir/annot_24hrsSS.txt $outdir/junction_stats_24hrs.bed $fasta $gtf &

regtools junctions annotate -S -o $outdir/annot_72hrsSS.txt $outdir/junction_stats_72hrs.bed $fasta $gtf ;

echo "Junction Annotations Complete"

