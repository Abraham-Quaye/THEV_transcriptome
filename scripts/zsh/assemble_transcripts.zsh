#!/usr/bin/env zsh


mapdir=results/hisat2
GTF=raw_files/annotations/thev_from_NCBI.gtf
assembled=results/stringtie


# constructing transcripts and splice variants with StringTie
echo "Assembling Sample Replicates..."

stringtie -p 10 $mapdir/thev_subset_72hrsS1.bam -t -a 1 -A $assembled/t72S1.tab -G $GTF -l 72S1 -o $assembled/thev_72hrsS1.gtf &
stringtie -p 10 $mapdir/thev_subset_72hrsS2.bam -t -a 1 -A $assembled/t72S2.tab -G $GTF -l 72S2 -o $assembled/thev_72hrsS2.gtf &&
stringtie -p 10 $mapdir/thev_subset_72hrsS3.bam -t -a 1 -g 1 -j 0.0001 -A $assembled/t72S3.tab -G $GTF -l 72S3 -o $assembled/thev_72hrsS3.gtf

stringtie -p 10 $mapdir/thev_subset_24hrsS1.bam -t -a 1 -A $assembled/t24S1.tab -G $GTF -l 24S1 -o $assembled/thev_24hrsS1.gtf &
stringtie -p 10 $mapdir/thev_subset_24hrsS2.bam -t -a 1 -A $assembled/t24S2.tab -G $GTF -l 24S2 -o $assembled/thev_24hrsS2.gtf &&
stringtie -p 10 $mapdir/thev_subset_24hrsS3.bam -t -a 1 -A $assembled/t24S3.tab -G $GTF -l 24S3 -o $assembled/thev_24hrsS3.gtf

stringtie -p 10 $mapdir/thev_subset_12hrsS1.bam -t -a 5 -A $assembled/t12S1.tab -G $GTF -l 12S1 -o $assembled/thev_12hrsS1.gtf &
stringtie -p 10 $mapdir/thev_subset_12hrsS3.bam -t -a 5 -A $assembled/t12S3.tab -G $GTF -l 12S3 -o $assembled/thev_12hrsS3.gtf &&

stringtie -p 10 $mapdir/thev_subset_4hrsS1.bam -t -A $assembled/t4S1.tab -G $GTF -l 4S1 -o $assembled/thev_4hrsS1.gtf &
stringtie -p 10 $mapdir/thev_subset_4hrsS2.bam -t -A $assembled/t4S2.tab -G $GTF -l 4S2 -o $assembled/thev_4hrsS2.gtf &
stringtie -p 10 $mapdir/thev_subset_4hrsS3.bam -t -A $assembled/t4S3.tab -G $GTF -l 4S3 -o $assembled/thev_4hrsS3.gtf &&

echo "All transcript assemblies complete!"

echo "Remove empty file: 4hrsS3"
rm $assembled/thev_4hrsS3.gtf $assembled/t4S3.tab

# echo "Create gtf merge list"
# ls $assembled/thev_4*.gtf > $assembled/merge_4hr_gtfs.txt && echo "4hr merge successful"
# ls $assembled/thev_12*.gtf > $assembled/merge_12hr_gtfs.txt && echo "12hr merge successful"
# ls $assembled/thev_24*.gtf > $assembled/merge_24hr_gtfs.txt && echo "24hr merge successful"
# ls $assembled/thev_72*.gtf > $assembled/merge_72hr_gtfs.txt && echo "72hr merge successful"


# # merge function 
# echo "Merging files..."

# stringtie --merge -i -t -j 0.0001 -f 0.0 -p 8 -o $assembled/final_4hr_trxpts.gtf -G $GTF $assembled/merge_4hr_gtfs.txt
# stringtie --merge -i -t -j 0.0001 -f 0.0 -p 8 -o $assembled/final_12hr_trxpts.gtf -G $GTF $assembled/merge_12hr_gtfs.txt
# stringtie --merge -i -t -j 0.0001 -f 0.0 -p 8 -o $assembled/final_24hr_trxpts.gtf -G $GTF $assembled/merge_24hr_gtfs.txt
# stringtie --merge -i -t -j 0.0001 -f 0.0 -p 8 -o $assembled/final_72hr_trxpts.gtf -G $GTF $assembled/merge_72hr_gtfs.txt

# echo "Stringtie script complete"
