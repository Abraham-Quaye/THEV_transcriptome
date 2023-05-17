#!/usr/bin/env zsh


mapdir=results/hisat2/bulk/uninfected
seqidx=raw_files/thevgenome_index/thev_tran
reads=trimmedReads/uninfected_reads

# map 72 hrs sample
(echo "mapping 72 hrs..."
hisat2 -p 10 --dta -x $seqidx -1 $reads/LCS9132_U_72hrsN1_Clean_Data1.fq.gz,$reads/LCS9132_U_72hrsN2_Clean_Data1.fq.gz -2 $reads/LCS9132_U_72hrsN1_Clean_Data2.fq.gz,$reads/LCS9132_U_72hrsN2_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/sortedTHEV_72hrsNeg.bam)

# map 24 hrs sample
(echo "mapping 24hrs..." 
hisat2 -p 10 --dta -x $seqidx -1 $reads/LCS9132_U_24hrsN1_Clean_Data1.fq.gz,$reads/LCS9132_U_24hrsN2_Clean_Data1.fq.gz -2 $reads/LCS9132_U_24hrsN1_Clean_Data2.fq.gz,$reads/LCS9132_U_24hrsN2_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/sortedTHEV_24hrsNeg.bam)

# map 12 hrs sample
(echo "mapping 12hrs..." 
hisat2 -p 10 --dta -x $seqidx -1 $reads/LCS9132_U_12hrsN1_Clean_Data1.fq.gz,$reads/LCS9132_U_12hrsN2_Clean_Data1.fq.gz -2 $reads/LCS9132_U_12hrsN1_Clean_Data2.fq.gz,$reads/LCS9132_U_12hrsN2_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/sortedTHEV_12hrsNeg.bam)

 # map 4 hrs sample
(echo "mapping 4hrs..." 
hisat2 -p 10 --dta -x $seqidx -1 $reads/LCS9132_U_4hrsN1_Clean_Data1.fq.gz,$reads/LCS9132_U_4hrsN2_Clean_Data1.fq.gz -2 $reads/LCS9132_U_4hrsN1_Clean_Data2.fq.gz,$reads/LCS9132_U_4hrsN2_Clean_Data2.fq.gz | samtools sort -@ 10 -o $mapdir/sortedTHEV_4hrsNeg.bam)

echo "Completed mapping at: $(date)"
