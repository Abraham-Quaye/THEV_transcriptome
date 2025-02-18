#!/usr/bin/env zsh

seqidx=raw_files/host_virus_genome_index/mix_tran
fordata=trimmedReads/forwardTrims
revdata=trimmedReads/reverseTrims
mapdir=results/hisat2


####################################### SECTION I ###############################################
# MAPPING ALL READS TO THEV GENOME WITH HISAT2

# map 72 hrs sample rep1
echo "Mapping 72hrs Replicate 1 ..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_72hrsS1_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_72hrsS1_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/thev_sorted_72hrsS1.bam

# map 72 hrs sample rep2
echo "Mapping 72hrs Replicate 2 ..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_72hrsS2_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_72hrsS2_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/thev_sorted_72hrsS2.bam

# map 72 hrs sample rep3
echo "Mapping 72hrs Replicate 3 ..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_72hrsS3_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_72hrsS3_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/thev_sorted_72hrsS3.bam

# map 24 hrs sample rep1
echo "Mapping 24hrs Replicate 1 ..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_24hrsS1_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_24hrsS1_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/thev_sorted_24hrsS1.bam

# map 24 hrs sample rep2
echo "Mapping 24hrs Replicate 2 ..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_24hrsS2_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_24hrsS2_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/thev_sorted_24hrsS2.bam

# map 24 hrs sample rep3
echo "Mapping 24hrs Replicate 3 ..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_24hrsS3_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_24hrsS3_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/thev_sorted_24hrsS3.bam

# map 12 hrs sample rep1
echo "Mapping 12hrs Replicate 1 ..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_12hrsS1_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_12hrsS1_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/thev_sorted_12hrsS1.bam

# map 12 hrs sample rep3
echo "Mapping 12hrs Replicate 3 ..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_12hrsS3_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_12hrsS3_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/thev_sorted_12hrsS3.bam

# map 4 hrs sample rep1
echo "Mapping 4hrs Replicate 1 ..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_4hrsS1_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_4hrsS1_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/thev_sorted_4hrsS1.bam

# map 4 hrs sample rep2
echo "Mapping 4hrs Replicate 2 ..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_4hrsS2_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_4hrsS2_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/thev_sorted_4hrsS2.bam

# map 4 hrs sample rep3
echo "Mapping 4hrs Replicate 3 ..."
hisat2 -p 10 --dta -x $seqidx -1 $fordata/LCS9132_I_4hrsS3_Clean_Data1_val_1.fq.gz -2 $revdata/LCS9132_I_4hrsS3_Clean_Data2_val_2.fq.gz | samtools sort -@ 10 -o $mapdir/thev_sorted_4hrsS3.bam

echo "All mapping complete!!!"
echo "Mapping Completed !!!"
