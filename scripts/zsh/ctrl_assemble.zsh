#!/usr/bin/env zsh


mapdir=results/hisat2/map_mock
GTF=raw_files/annotations/thev_from_NCBI.gtf
assembled=results/stringtie/mock_stringtie


############################################# SECTION V #############################################
# constructing transcripts and splice variants with StringTie
echo "Assembling ctrl Replicates..."

stringtie -p 10 $mapdir/ctrl_sorted_72N1.bam -t -G $GTF -l 72S1 -o $assembled/ctrl_72N1.gtf &
stringtie -p 10 $mapdir/ctrl_sorted_72N2.bam -t -G $GTF -l 72S2 -o $assembled/ctrl_72N2.gtf &

stringtie -p 10 $mapdir/ctrl_sorted_24N1.bam -t -G $GTF -l 24S1 -o $assembled/ctrl_24N1.gtf &
stringtie -p 10 $mapdir/ctrl_sorted_24N2.bam -t -G $GTF -l 24S2 -o $assembled/ctrl_24N2.gtf &

stringtie -p 10 $mapdir/ctrl_sorted_12N1.bam -t -G $GTF -l 12S1 -o $assembled/ctrl_12N1.gtf &
stringtie -p 10 $mapdir/ctrl_sorted_12N2.bam -t -G $GTF -l 12S2 -o $assembled/ctrl_12N2.gtf &&

stringtie -p 10 $mapdir/ctrl_sorted_4N1.bam -t -G $GTF -l 4S1 -o $assembled/ctrl_4N1.gtf &
stringtie -p 10 $mapdir/ctrl_sorted_4N2.bam -t -G $GTF -l 4S2 -o $assembled/ctrl_4N2.gtf ;

echo "All ctrl transcript assemblies complete!"

# merge function 
### stringtie --merge -i -f -p 8 -o merge_imp_all.gtf ./merge_imp_all.txt

