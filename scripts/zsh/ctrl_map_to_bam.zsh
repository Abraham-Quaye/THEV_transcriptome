#!/usr/bin/env zsh

seqidx=raw_files/thevgenome_index/thev_tran
data=trimmedReads/uninfected_reads
mapdir=results/hisat2/map_mock

####################################### SECTION I ###############################################
# MAPPING ALL READS TO THEV GENOME WITH HISAT2

# map 4 hrs sample rep1
(echo "Mapping ctrl 4hrs Replicate 1 ..."
hisat2 -p 10 --dta -x $seqidx -1 $data/LCS9132_U_4hrsN1_Clean_Data1.fq.gz -2 $data/LCS9132_U_4hrsN1_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/ctrl_sorted_4N1.bam) &

# map 4 hrs sample rep2
(echo "Mapping ctrl 4hrs Replicate 2 ..."
hisat2 -p 10 --dta -x $seqidx -1 $data/LCS9132_U_4hrsN2_Clean_Data1.fq.gz -2 $data/LCS9132_U_4hrsN2_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/ctrl_sorted_4N2.bam) ;

# map 72 hrs sample rep1
(echo "Mapping ctrl 72hrs Replicate 1 ..."
hisat2 -p 10 --dta -x $seqidx -1 $data/LCS9132_U_72hrsN1_Clean_Data1.fq.gz -2 $data/LCS9132_U_72hrsN1_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/ctrl_sorted_72N1.bam) &

# map 72 hrs sample rep2
(echo "Mapping ctrl 72hrs Replicate 2 ..."
hisat2 -p 10 --dta -x $seqidx -1 $data/LCS9132_U_72hrsN2_Clean_Data1.fq.gz -2 $data/LCS9132_U_72hrsN2_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/ctrl_sorted_72N2.bam) ;

# map 24 hrs sample rep1
(echo "Mapping ctrl 24hrs Replicate 1 ..."
hisat2 -p 10 --dta -x $seqidx -1 $data/LCS9132_U_24hrsN1_Clean_Data1.fq.gz -2 $data/LCS9132_U_24hrsN1_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/ctrl_sorted_24N1.bam) &

# map 24 hrs sample rep2
(echo "Mapping ctrl 24hrs Replicate 2 ..."
hisat2 -p 10 --dta -x $seqidx -1 $data/LCS9132_U_24hrsN2_Clean_Data1.fq.gz -2 $data/LCS9132_U_24hrsN2_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/ctrl_sorted_24N2.bam) ;

# map 12 hrs sample rep1
(echo "Mapping ctrl 12hrs Replicate 1 ..."
hisat2 -p 10 --dta -x $seqidx -1 $data/LCS9132_U_12hrsN1_Clean_Data1.fq.gz -2 $data/LCS9132_U_12hrsN1_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/ctrl_sorted_12N1.bam) &

# map 12 hrs sample rep3
(echo "Mapping ctrl 12hrs Replicate 2 ..."
hisat2 -p 10 --dta -x $seqidx -1 $data/LCS9132_U_12hrsN2_Clean_Data1.fq.gz -2 $data/LCS9132_U_12hrsN2_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/ctrl_sorted_12N2.bam) ;

echo "All ctrl mapping complete!!!"
