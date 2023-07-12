#!/usr/bin/env zsh


mapdir=results/hisat2/bulk
seqidx=raw_files/host_virus_genome_index/mix_tran
fordata=trimmedReads/forwardTrims
revdata=trimmedReads/reverseTrims

# map 72 hrs sample
(echo "mapping 72 hrs..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_72hrsS1_Clean_Data1_val_1.fq.gz,$fordata/LCS9132_I_72hrsS2_Clean_Data1_val_1.fq.gz,$fordata/LCS9132_I_72hrsS3_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_72hrsS1_Clean_Data2_val_2.fq.gz,$revdata/LCS9132_I_72hrsS2_Clean_Data2_val_2.fq.gz,$revdata/LCS9132_I_72hrsS3_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/sortedTHEV_72hrsSamples.bam) ;

# map 24 hrs sample
(echo "mapping 24hrs..." 
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_24hrsS1_Clean_Data1_val_1.fq.gz,$fordata/LCS9132_I_24hrsS2_Clean_Data1_val_1.fq.gz,$fordata/LCS9132_I_24hrsS3_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_24hrsS1_Clean_Data2_val_2.fq.gz,$revdata/LCS9132_I_24hrsS2_Clean_Data2_val_2.fq.gz,$revdata/LCS9132_I_24hrsS3_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/sortedTHEV_24hrsSamples.bam) ;

# map 12 hrs sample
(echo "mapping 12hrs..." 
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_12hrsS1_Clean_Data1_val_1.fq.gz,$fordata/LCS9132_I_12hrsS3_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_12hrsS1_Clean_Data2_val_2.fq.gz,$revdata/LCS9132_I_12hrsS3_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/sortedTHEV_12hrsSamples.bam) ;

 # map 4 hrs sample
(echo "mapping 4hrs..." 
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_4hrsS1_Clean_Data1_val_1.fq.gz,$fordata/LCS9132_I_4hrsS2_Clean_Data1_val_1.fq.gz,$fordata/LCS9132_I_4hrsS3_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_4hrsS1_Clean_Data2_val_2.fq.gz,$revdata/LCS9132_I_4hrsS2_Clean_Data2_val_2.fq.gz,$revdata/LCS9132_I_4hrsS3_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/sortedTHEV_4hrsSamples.bam) ;

echo "Completed mapping at: $(date)"
