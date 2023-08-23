
################# CONVERT .BED TO GFF AND MOVE AGAT LOGFILE ###########
rule make_gff:
    input:
        script = "scripts/zsh/make_gff.zsh",
        r_script = "scripts/r/trim_gff.R",
        bed = "raw_files/annotations/THEVannotated_genesOnly.bed"
    output:
        "raw_files/annotations/thev_predicted_genes.gff"
    shell:
        """
        {input.script}
        Rscript scripts/r/trim_gff.R
        """

################# CONVERT .GFF3 TO GTF AND MOVE AGAT LOGFILE ###########
rule make_gtf:
    input:
        script = "scripts/zsh/make_gtf.zsh",
        thevgff = rules.make_gff.output,
        hostgff = "raw_files/annotations/turkey_genome.gff"
    output:
        "raw_files/annotations/thev_predicted_genes.gtf",
        "raw_files/annotations/thev_predicted_genes.agat.log",
        "raw_files/annotations/turkey_genome.gtf",
        "raw_files/annotations/turkey_genome.agat.log"
    shell:
        "{input.script}"

################### EXTRACT SPLICE-SITES ####################
rule extract_splice_site:
    input:
        script = "scripts/zsh/extract_ss.zsh",
        thevgtf = "raw_files/annotations/thev_predicted_genes.gtf",
        hostgtf = "raw_files/annotations/turkey_genome.gtf"
    output:
        "raw_files/annotations/thev_predicted_genes.ss",
        "raw_files/annotations/turkey_genome.ss",
        "raw_files/annotations/host_thev.ss"
    shell:
        "{input.script}"

################### EXTRACT EXONS ####################
rule extract_exons:
    input:
        bash_script = "scripts/zsh/extract_exons.zsh",
        r_script = "scripts/r/extract_thev_exons.R",
        thevgtf = "raw_files/annotations/thev_predicted_genes.gtf",
        hostgtf = "raw_files/annotations/turkey_genome.gtf"
    output:
        "raw_files/annotations/thev_predicted_genes.exons",
        "raw_files/annotations/turkey_genome.exons",
        "raw_files/annotations/host_thev.exons"
    shell:
        """
        {input.r_script}
        {input.bash_script}
        """

#################### BUILD THEV GENOMIC INDEX FOR MAPPING WITH HISAT2 ######
rule build_genome_index:
    input:
        script = "scripts/zsh/build_genome_index.zsh",
        ss = "raw_files/annotations/host_thev.ss",
        exon = "raw_files/annotations/host_thev.exons",
        mix_genome = "raw_files/genome_file/mix_genomes.fa"
    output:
        expand("raw_files/host_virus_genome_index/mix_tran.{n}.ht2", \
        n = range(1,9))
    shell:
        "{input.script}"

#################### MAP READS TO THEV GENOME WITH HISAT2 #############
rule map_sort_to_bam:
    input:
        script = "scripts/zsh/map_sort_to_bam.zsh",
        seqidx = rules.build_genome_index.output,
        fordata = expand("trimmedReads/forwardTrims/LCS9132_I_{tp}hrsS{rep}_Clean_Data1_val_1.fq.gz", \
        tp = [72, 24, 4], rep = [1, 2, 3]),
        for12 = expand("trimmedReads/forwardTrims/LCS9132_I_12hrsS{rep}_Clean_Data1_val_1.fq.gz", \ 
        rep = [1, 3] ),
        revdata = expand("trimmedReads/reverseTrims/LCS9132_I_{tp}hrsS{rep}_Clean_Data2_val_2.fq.gz", \
        tp = [72, 24, 4], rep = [1, 2, 3]),
        rev12 = expand("trimmedReads/reverseTrims/LCS9132_I_12hrsS{rep}_Clean_Data2_val_2.fq.gz", \
        rep = [1, 3])
    output:
        expand("results/hisat2/thev_sorted_{time}hrsS{rep}.bam", \
        time = ["72", "24", "4"], rep = ["1", "2", "3"]),
        expand("results/hisat2/thev_sorted_12hrsS{rep}.bam", \
        rep = ["1", "3"])
    shell:
        "{input.script}"

#################### INDEX ALL SORTED .BAM FILES ##################
rule index_bam_files:
    input:
        script = "scripts/zsh/index.zsh",
        bam = rules.map_sort_to_bam.output
    output:
        expand("results/hisat2/thev_sorted_{time}hrsS{rep}.bam.bai", \
        time = ["72", "24", "4"], rep = ["1", "2", "3"]),
        expand("results/hisat2/thev_sorted_12hrsS{rep}.bam.bai", \
        rep = ["1", "3"])
    shell:
        "{input.script}"

#################### SUBSET SORTED BAM FILES FOR THEV ##################
rule filter_thev:
    input:
        script = "scripts/zsh/thev_subset.zsh",
        bam = rules.map_sort_to_bam.output,
        idx = rules.index_bam_files.output
    output:
        expand("results/hisat2/thev_subset_{time}hrsS{rep}.bam", \
        time = ["72", "24", "4"], rep = ["1", "2", "3"]),
        expand("results/hisat2/thev_subset_12hrsS{rep}.bam", \
        rep = ["1", "3"])
    shell:
        "{input.script}"

#################### INDEX ALL SUBSET BAM FILES ##################
rule index_subset_bam:
    input:
        script = "scripts/zsh/subset_index.zsh",
        bam = rules.filter_thev.output
    output:
        expand("results/hisat2/thev_subset_{time}hrsS{rep}.bam.bai", \
        time = ["72", "24", "4"], rep = ["1", "2", "3"]),
        expand("results/hisat2/thev_subset_12hrsS{rep}.bam.bai", \
        rep = ["1", "3"])
    shell:
        "{input.script}"

#################### CONSTRUCT TRANSCRIPTS WITH STRINGTIE ##############
rule make_transcripts:
    input:
        script = "scripts/zsh/assemble_transcripts.zsh",
        gtf = "raw_files/annotations/thev_predicted_genes.gtf",
        bam = rules.filter_thev.output
    output:
        expand("results/stringtie/thev_{time}hrsS{rep}.gtf", \
        time = [24, 72], rep = [1, 2, 3]),
        expand("results/stringtie/thev_12hrsS{rep}.gtf", \
        rep = [1, 3]),
        expand("results/stringtie/thev_4hrsS{rep}.gtf", \
        rep = [1, 2]),
        expand("results/stringtie/t{time}S{rep}.tab", \
        time = [24, 72], rep = [1, 2, 3]),
        expand("results/stringtie/t12S{rep}.tab", \
        rep = [1, 3]),
        expand("results/stringtie/t4S{rep}.tab", \
        rep = [1, 2])
    shell:
        "{input.script}"

#################### MERGE ALL GTF FILES #####################
rule merge_gtfs:
    input:
        expand("results/stringtie/thev_{time}hrsS{rep}.gtf", \
        time = [24, 72], rep = [1, 2, 3]),
        expand("results/stringtie/thev_12hrsS{rep}.gtf", \
        rep = [1, 3]),
        expand("results/stringtie/thev_4hrsS{rep}.gtf", \
        rep = [1, 2])
    output:
        "results/stringtie/all_merged.gtf"
    shell:
        "cat {input} > {output}"

##################### FILTER FOR REAL TRANSCRIPTS (REMOVE PREDICTED ORFs) ######
rule split_timepoint_transcripts:
    input:
        r_script = "scripts/r/filter_real_transcripts.R",
        all_gtf = rules.merge_gtfs.output
    output:
        "results/stringtie/all_real_transcripts_merged.gtf",
        expand("results/stringtie/transcripts_merged_{tp}hrs.gtf", \
        tp = [4, 12, 24, 72])
    shell:
        "{input.r_script}"

############## GFFCOMPARE TO GENERATE FINAL TRANSCRIPTOME ##############
rule make_unredundant_transcriptome:
    input:
        gtfs = expand("results/stringtie/transcripts_merged_{tp}hrs.gtf", \
        tp = [4, 12, 24, 72]),
        script = "scripts/zsh/gffcompare.zsh",
        orfs = "raw_files/annotations/thev_predicted_genes.gtf"
    output:
        expand("results/gffcompare/gffcomp_alltimes.{metric}", \
        metric = ["combined.gtf", "loci", "stats", "tracking"])
    shell:
        """
        {input.script}
        mv results/stringtie/gffcomp_alltimes.* results/gffcompare/
        """

############## ADD REGION ATTRIBUTE TO GFFCOMPARE GTF FILE ##############
rule mod_final_trxptome:
    input:
        gtf = "results/gffcompare/gffcomp_alltimes.combined.gtf",
        rscript = "scripts/r/add_region_to_trxptome.R"
    output:
        "results/gffcompare/updated_alltimes.combined.gtf"
    shell:
        "{input.rscript}"

################### MAKE FULL TRANSCRIPTOME WITH UNREDUNDANT GTF FILE FROM GFFCOMPARE ############
rule make_full_splice_map:
    input:
        spliced_gtf = "results/gffcompare/gffcomp_alltimes.combined.gtf",
        orf_gtf = rules.merge_gtfs.output,
        rscript = "scripts/r/thev_splicing_fullMap.R"
    output:
        "results/r/figures/thev_spliced_map.png"
    shell:
        "{input.rscript}"

################### MAKE TRANSCRIPTOME MAP PER TIMEPOINT OF INFECTION #####################
rule make_timepoint_splice_map:
    input:
        rscript = "scripts/r/plot_thev_timepoint_splicing.R",
        gtfs = expand("results/stringtie/transcripts_merged_{tp}hrs.gtf", \
        tp = [4, 12, 24, 72])
    output:
        "results/r/figures/thev_patched_timepoints_spliced_map.png"
    shell:
        "{input.rscript}"

############## ESTIMATE TRANSCRIPT ABUNDANCES##############
rule est_abundances:
    input:
        script = "scripts/zsh/est_trxpt_abund.zsh",
        bams = rules.filter_thev.output,
        ref = "results/stringtie/all_real_transcripts_merged.gtf",
    output:
        expand("results/ballgown/abund_{s}/abund_{s}.gtf", \
        s = ["72hrsS1", "72hrsS2", "72hrsS3", "24hrsS1", "24hrsS2", "24hrsS3", "12hrsS1", "12hrsS3", "4hrsS1", "4hrsS2"])
    shell:
        "{input.script}"

#################### COUNT ALL SPLICE JUNCTIONS ####################
rule count_junctions:
    input:
        bams = rules.filter_thev.output,
        index = rules.index_subset_bam.output,
        script = "scripts/zsh/splice_site_stats.zsh"
    output:
        expand("results/hisat2/junction_stats_{tp}S{rep}.bed", tp = [4, 24, 72], rep = [1, 2, 3]),
        expand("results/hisat2/junction_stats_12S{rep}.bed", rep = [1, 3])
    shell:
        "{input.script}"

#################### BULK MAP READS ##########
rule bulk_map_sort_to_bam:
    input:
        script = "scripts/zsh/bulk_map_sort_to_bam.zsh",
        seqidx = rules.build_genome_index.output,
        fordata = expand("trimmedReads/forwardTrims/LCS9132_I_{tp}hrsS{rep}_Clean_Data1_val_1.fq.gz", \
        tp = [72, 24, 4], rep = [1, 2, 3]),
        for12 = expand("trimmedReads/forwardTrims/LCS9132_I_12hrsS{rep}_Clean_Data1_val_1.fq.gz", \ 
        rep = [1, 3]),
        revdata = expand("trimmedReads/reverseTrims/LCS9132_I_{tp}hrsS{rep}_Clean_Data2_val_2.fq.gz", \
        tp = [72, 24, 4], rep = [1, 2, 3]),
        rev12 = expand("trimmedReads/reverseTrims/LCS9132_I_12hrsS{rep}_Clean_Data2_val_2.fq.gz", \
        rep = [1, 3])
    output:
        expand("results/hisat2/bulk/sortedTHEV_{time}hrsSamples.bam", \
        time = [4, 12, 24, 72])
    shell:
        "{input.script}"

#################### BULK INDEX ############
rule bulk_index:
    input:
        script = "scripts/zsh/bulk_index.zsh",
        bam = rules.bulk_map_sort_to_bam.output
    output:
        expand("results/hisat2/bulk/sortedTHEV_{time}hrsSamples.bam.bai", \
        time = [4, 12, 24, 72])
    shell:
        "{input.script}"

#################### HOST AND VIRUS BULK COVERAGE #########
rule total_bulk_coverage:
    input:
        script = "scripts/zsh/total_bulk_coverage.zsh",
        bam = rules.bulk_map_sort_to_bam.output,
        idx = rules.bulk_index.output
    output:
        expand("results/hisat2/coverage/host_thev_cov{t}.txt", \
        t = [4, 12, 24, 72])    
    shell:
        "{input.script}"

#################### BULK SUBSET THEV #########
rule filter_bulk_thev:
    input:
        script = "scripts/zsh/subset_bulk_thev.zsh",
        bams = rules.bulk_map_sort_to_bam.output,
        idx = rules.bulk_index.output
    output:
        expand("results/hisat2/bulk/subsetTHEV_{time}hrs.bam", \
        time = [4, 12, 24, 72])
    shell:
        "{input.script}"

#################### BULK SUBSET INDEX ############
rule bulk_subset_index:
    input:
        script = "scripts/zsh/subset_bulk_idx.zsh",
        bams = rules.filter_bulk_thev.output
    output:
        expand("results/hisat2/bulk/subsetTHEV_{time}hrs.bam.bai", \
        time = [4, 12, 24, 72])
    shell:
        "{input.script}"

#################### BULK DEPTH ############
rule bulk_depth:
    input:
        script = "scripts/zsh/bulk_depth.zsh",
        bam = rules.filter_bulk_thev.output,
        idx = rules.bulk_subset_index.output
    output:
        expand("results/hisat2/coverage/thev_{time}hrsdepth.txt", \
        time = [4, 12, 24, 72])
    shell:
        "{input.script}"

#################### BULK COUNT JUNCTIONS ############
rule bulk_count_junctions:
    input:
        bams = rules.filter_bulk_thev.output,
        index = rules.bulk_subset_index.output,
        script = "scripts/zsh/bulk_splice_stats.zsh"
    output:
        expand("results/hisat2/bulk/junction_stats_{tp}hrs.bed", tp = [4, 12, 24, 72])
    shell:
        "{input.script}"

#################### BULK ANNOTATE JUNCTIONS ############
rule bulk_annotate_junctions:
    input:
        script = "scripts/zsh/annotate_junctions.zsh",
        bedfiles = rules.bulk_count_junctions.output,
        fasta = "raw_files/genome_file/AY849321.1.fa",
        gtf= rules.mod_final_trxptome.output
    output:
       expand("results/hisat2/bulk/annot_{tp}hrsSS.txt", tp = [4, 12, 24, 72])
    shell:
        "{input.script}"

#################### BULK COUNTING READS ############
rule count_total_reads:
    input:
        script = "scripts/zsh/count_total_reads.zsh",
        bam = rules.bulk_map_sort_to_bam.output
    output:
        "results/hisat2/coverage/bulk_counts.txt"
    shell:
        "{input.script}"

#################### MAKE FIGURES FOR DEPTH/COVERAGE ############
rule make_coverage_figures:
    input:
        r_script1 = "scripts/r/thev_cov_depth.R",
        r_script2 = "scripts/r/thev_genomic_map.R",
        bedfile = "raw_files/annotations/THEVannotated_genesOnly.bed",
        depth = rules.bulk_depth.output,
        coverage = rules.total_bulk_coverage.output
    output:
        # expand("results/r/figures/depth_{time}hrs.png", \
        # time = [4, 12, 24, 72]),
        expand("results/r/figures/{plotkind}_alltimes.png", \
        plotkind = ["patch", "overlay", "correlate"])
    shell:
       "{input.r_script1}"

################## PLOT JUNCTION ABUNDANCE FIGURES ################
rule junc_abund_plots:
    input:
        trxptome = rules.mod_final_trxptome.output,
        annot_juncs = rules.bulk_annotate_junctions.output,
        rscript = "scripts/r/splice_site_analyses.R"
    output:
        "results/r/figures/junc_abundances.png"
    shell:
        "{input.rscript}"

#################### PLOT THEV GENOMIC MAP ############
rule make_orf_map:
    input:
        r_script = "scripts/r/thev_orf_map.R",
        bedfile = "raw_files/annotations/THEVannotated_genesOnly.bed",
    output:
        "results/r/figures/thev_orf_map.png"
    shell:
       "{input.r_script}"

################# MAP UNINFECTED READS ############################
# rule uninfected_map:
#     input:
#         script = "scripts/zsh/ctrl_map_to_bam.zsh",
#         seqidx = expand("raw_files/thevgenome_index/thev_tran.{n}.ht2", \
#         n = [1, 2, 3, 4, 5, 6, 7, 8]),
#         reads = expand("trimmedReads/uninfected_reads/LCS9132_U_{tp}hrsN{rep}_Clean_Data{strand}.fq.gz", \
#         tp = [72, 24, 12, 4], rep = [1, 2], strand = [1, 2])
#     output:
#         expand("results/hisat2/map_mock/ctrl_sorted_{time}N{rep}.bam", \
#         time = [72, 24, 12, 4], rep = [1, 2])
#     shell:
#         "{input.script}"

# ################# INDEX UNINFECTED READS ############################
# rule uninfected_single_index:
#     input:
#         bam = rules.uninfected_map.output,
#         script = "scripts/zsh/ctrl_single_index.zsh"
#     output:
#         expand("results/hisat2/map_mock/ctrl_sorted_{time}N{rep}.bam.bai", \
#         time = ["72", "24", "12", "4"], rep = ["1", "2"])
#     shell:
#         "{input.script}"

# ################# ASSEMBLE UNINFECTED TRANSCRITPS ############################
# rule uninfected_assemble:
#     input:
#         bam = rules.uninfected_map.output,
#         script = "scripts/zsh/ctrl_assemble.zsh",
#         gtf = "raw_files/annotations/thev_from_NCBI.gtf"
#     output:
#         expand("results/stringtie/mock_stringtie/ctrl_{time}N{rep}.gtf", \
#         time = [4, 12, 24, 72], rep = [1, 2])
#     shell:
#         "{input.script}"

# ################# BULK MAP UNINFECTED READS ############################
# rule uninfected_map_to_bam:
#     input:
#         script = "scripts/zsh/ctrl_bulk_map_sort_to_bam.zsh",
#         seqidx = expand("raw_files/thevgenome_index/thev_tran.{n}.ht2", \
#         n = [1, 2, 3, 4, 5, 6, 7, 8]),
#         reads = expand("trimmedReads/uninfected_reads/LCS9132_U_{tp}hrsN{rep}_Clean_Data{strand}.fq.gz", \
#         tp = [72, 24, 12, 4], rep = [1, 2], strand = [1, 2])
#     output:
#         expand("results/hisat2/bulk/uninfected/sortedTHEV_{time}hrsNeg.bam", \
#         time = [4, 12, 24, 72])
#     shell:
#         "{input.script}"

# ############### INDEX UNINFECTED ##################
# rule uninfected_index:
#     input:
#         rules.uninfected_map_to_bam.output,
#         script = "scripts/zsh/ctrl_index.zsh"
#     output:
#         expand("results/hisat2/bulk/uninfected/sortedTHEV_{time}hrsNeg.bam.bai", \
#         time = [4, 12, 24, 72])
#     shell:
#         "{input.script}"

# ############### COVERAGE UNINFECTED ##################
# rule uninfected_coverage:
#     input:
#         rules.uninfected_map_to_bam.output,
#         script = "scripts/zsh/ctrl_coverage.zsh"
#     output:
#         "results/hisat2/coverage/ctrl_coverage.tsv"
#     shell:
#         "{input.script}"

# ############### DEPTH UNINFECTED ##################
# rule uninfected_depth:
#     input:
#         rules.uninfected_map_to_bam.output,
#         script = "scripts/zsh/ctrl_depth.zsh"
#     output:
#         expand("results/hisat2/coverage/ctrl_{time}hrsdepth.txt", \
#         time = [4, 12, 24, 72])
#     shell:
#         "{input.script}"

# ############### FIGURES FOR UNINFECTED ##################
# rule uninfected_coverage_figures:
#     input:
#         r_script1 = "scripts/r/ctrl_cov_depth.R",
#         r_script2 = "scripts/r/thev_genomic_map.R",
#         bedfile = "raw_files/annotations/THEVannotated_genesOnly.bed",
#         depth = expand("results/hisat2/coverage/ctrl_{time}hrsdepth.txt", \
#         time = [4, 12, 24, 72]),
#         coverage = "results/hisat2/coverage/ctrl_coverage.tsv"
#     output:
#         expand("results/r/figures/ctrl_depth_{time}hrs.pdf", \
#         time = [4, 12, 24, 72]),
#         expand("results/r/figures/ctrl_{plotkind}_alltimes.pdf", \
#         plotkind = ["patch", "overlay", "correlate"])
#     shell:
#        "{input.r_script1}"

#################### MAKE THEV GROWTH CURVE ############
rule make_growth_curve:
    input:
        r_script = "scripts/r/thev_growthcurve.R",
        exp1 = "raw_files/wetlab_data/thev_growthcurve2021.xls",
        exp2 = "raw_files/wetlab_data/thev_growthcurve04_2023.xls"
    output:
        "results/r/figures/thev_growth_curve.png"
    shell:
       "{input.r_script}"

############### WRITE MANUSCRIPT FOR PUBLICATION ##################
rule write_manuscript:
    input:
        "manuscript_thev_transcriptome.Rmd",
        "asm.csl",
        "transcriptome_refs.bib",
        "scripts/r/bam_file_analysis.R",
        rules.bulk_annotate_junctions.output,
        rules.total_bulk_coverage.output,
        rules.count_total_reads.output,
        rules.make_coverage_figures.output,
        rules.make_growth_curve.output,
        rules.make_orf_map.output,
        rules.make_full_splice_map.output,
        rules.make_timepoint_splice_map.output,
        rules.junc_abund_plots.output
    output:
        "manuscript_thev_transcriptome.html",
        # "manuscript_thev_transcriptome.docx"
    shell:
        """
        R -e "library(rmarkdown);render('manuscript_thev_transcriptome.Rmd', output_format = 'all')"
        """

rule write_supplementary:
    input:
        "supplementary_thev_trxptome.Rmd",
        "scripts/r/bam_file_analysis.R",
    output:
        "manuscript_thev_transcriptome.pdf"
    shell:
        """
        R -e "library(rmarkdown);render('supplementary_thev_trxptome.Rmd')"
        """
############# RUN ENTIRE SCRIPT RULE ##############
rule run_pipeline:
    input:
        rules.est_abundances.output,
        rules.write_manuscript.output,
        rules.write_supplementary.output
        