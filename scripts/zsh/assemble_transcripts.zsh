#!/usr/bin/env zsh


mapdir=results/hisat2
GTF=raw_files/annotations/thev_predicted_genes.gtf
assembled=results/stringtie


# constructing transcripts and splice variants with StringTie
echo "Assembling Sample Replicates..."

stringtie -p 10 $mapdir/thev_subset_72hrsS1.bam -t -a 1 -G $GTF -l 72S1 -o $assembled/thev_72hrsS1.gtf &
stringtie -p 10 $mapdir/thev_subset_72hrsS2.bam -t -a 1 -G $GTF -l 72S2 -o $assembled/thev_72hrsS2.gtf &&
stringtie -p 10 $mapdir/thev_subset_72hrsS3.bam -t -a 1 -g 1 -j 0.0001 -G $GTF -l 72S3 -o $assembled/thev_72hrsS3.gtf

stringtie -p 10 $mapdir/thev_subset_24hrsS1.bam -t -a 1 -f 0.02 -G $GTF -l 24S1 -o $assembled/thev_24hrsS1.gtf &
stringtie -p 10 $mapdir/thev_subset_24hrsS2.bam -t -a 1 -G $GTF -l 24S2 -o $assembled/thev_24hrsS2.gtf &&
stringtie -p 10 $mapdir/thev_subset_24hrsS3.bam -t -a 1 -G $GTF -l 24S3 -o $assembled/thev_24hrsS3.gtf

stringtie -p 10 $mapdir/thev_subset_12hrsS1.bam -t -a 5 -G $GTF -l 12S1 -o $assembled/thev_12hrsS1.gtf &
stringtie -p 10 $mapdir/thev_subset_12hrsS3.bam -t -a 5 -G $GTF -l 12S3 -o $assembled/thev_12hrsS3.gtf &&

stringtie -p 10 $mapdir/thev_subset_4hrsS1.bam -t -G $GTF -l 4S1 -o $assembled/thev_4hrsS1.gtf &
stringtie -p 10 $mapdir/thev_subset_4hrsS2.bam -t -G $GTF -l 4S2 -o $assembled/thev_4hrsS2.gtf &
stringtie -p 10 $mapdir/thev_subset_4hrsS3.bam -t -G $GTF -l 4S3 -o $assembled/thev_4hrsS3.gtf &&

echo "All transcript assemblies complete!"

echo "Remove empty file: 4hrsS3"
rm $assembled/thev_4hrsS3.gtf

