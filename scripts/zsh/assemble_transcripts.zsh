#!/usr/bin/env zsh


mapdir=results/hisat2
GFF=raw_files/annotations/thev_predicted_genes.gff
assembled=results/stringtie
############################################# SECTION V #############################################
# constructing transcripts and splice variants with StringTie
echo "SECTION V beginning..."
echo "Assembling Sample Replicates..."

stringtie -p 10 $mapdir/thev_sorted_4hrsS1.bam -t -G $GFF -l 4S1 -o $assembled/thev_4hrsS1.gtf &
stringtie -p 10 $mapdir/thev_sorted_4hrsS2.bam -t -G $GFF -l 4S2 -o $assembled/thev_4hrsS2.gtf &
stringtie -p 10 $mapdir/thev_sorted_4hrsS3.bam -t -G $GFF -l 4S3 -o $assembled/thev_4hrsS3.gtf &

stringtie -p 10 $mapdir/thev_sorted_12hrsS1.bam -t -G $GFF -l 12S1 -o $assembled/thev_12hrsS1.gtf &
stringtie -p 10 $mapdir/thev_sorted_12hrsS3.bam -t -G $GFF -l 12S3 -o $assembled/thev_12hrsS3.gtf &

stringtie -p 10 $mapdir/thev_sorted_24hrsS1.bam -t -G $GFF -l 24S1 -o $assembled/thev_24hrsS1.gtf &
stringtie -p 10 $mapdir/thev_sorted_24hrsS2.bam -t -G $GFF -l 24S2 -o $assembled/thev_24hrsS2.gtf &
stringtie -p 10 $mapdir/thev_sorted_24hrsS3.bam -t -G $GFF -l 24S3 -o $assembled/thev_24hrsS3.gtf &

stringtie -p 10 $mapdir/thev_sorted_72hrsS1.bam -t -G $GFF -l 72S1 -o $assembled/thev_72hrsS1.gtf &
stringtie -p 10 $mapdir/thev_sorted_72hrsS2.bam -t -G $GFF -l 72S2 -o $assembled/thev_72hrsS2.gtf &
stringtie -p 10 $mapdir/thev_sorted_72hrsS3.bam -t -G $GFF -l 72S3 -o $assembled/thev_72hrsS3.gtf ;

echo "All transcript assemblies complete!"

# merge function 
### stringtie --merge -i -f -p 8 -o merge_imp_all.gtf ./merge_imp_all.txt