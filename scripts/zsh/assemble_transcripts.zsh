#!/usr/bin/env zsh


mapdir=results/hisat2
GTF=raw_files/annotations/thev_from_NCBI.gtf
assembled=results/stringtie

rm $assembled/*.gtf && echo "old .gtf files removed" || echo "no .gtf files in directory"

############################################# SECTION V #############################################
# constructing transcripts and splice variants with StringTie
echo "Assembling Sample Replicates..."

stringtie -p 10 $mapdir/thev_sorted_72hrsS1.bam -t -A $assembled/t72S1.tab -G $GTF -l 72S1 -o $assembled/thev_72hrsS1.gtf &
stringtie -p 10 $mapdir/thev_sorted_72hrsS2.bam -t -A $assembled/t72S2.tab -G $GTF -l 72S2 -o $assembled/thev_72hrsS2.gtf &
stringtie -p 10 $mapdir/thev_sorted_72hrsS3.bam -t -A $assembled/t72S3.tab -G $GTF -l 72S3 -o $assembled/thev_72hrsS3.gtf &

stringtie -p 10 $mapdir/thev_sorted_24hrsS1.bam -t -A $assembled/t24S1.tab -G $GTF -l 24S1 -o $assembled/thev_24hrsS1.gtf &
stringtie -p 10 $mapdir/thev_sorted_24hrsS2.bam -t -A $assembled/t24S2.tab -G $GTF -l 24S2 -o $assembled/thev_24hrsS2.gtf &
stringtie -p 10 $mapdir/thev_sorted_24hrsS3.bam -t -A $assembled/t24S3.tab -G $GTF -l 24S3 -o $assembled/thev_24hrsS3.gtf &

stringtie -p 10 $mapdir/thev_sorted_12hrsS1.bam -t -A $assembled/t12S1.tab -G $GTF -l 12S1 -o $assembled/thev_12hrsS1.gtf &
stringtie -p 10 $mapdir/thev_sorted_12hrsS3.bam -t -A $assembled/t12S3.tab -G $GTF -l 12S3 -o $assembled/thev_12hrsS3.gtf &&

stringtie -p 10 $mapdir/thev_sorted_4hrsS1.bam -t -A $assembled/t4S1.tab -G $GTF -l 4S1 -o $assembled/thev_4hrsS1.gtf &
stringtie -p 10 $mapdir/thev_sorted_4hrsS2.bam -t -A $assembled/t4S2.tab -G $GTF -l 4S2 -o $assembled/thev_4hrsS2.gtf &
stringtie -p 10 $mapdir/thev_sorted_4hrsS3.bam -t -A $assembled/t4S3.tab -G $GTF -l 4S3 -o $assembled/thev_4hrsS3.gtf &&

echo "All transcript assemblies complete!"

# merge function 
### stringtie --merge -i -f -p 8 -o merge_imp_all.gtf ./merge_imp_all.txt
