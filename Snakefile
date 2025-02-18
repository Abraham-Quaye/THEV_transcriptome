
################# CONVERT .BED TO GFF AND MOVE AGAT LOGFILE ###########
rule make_gff:
    input:
        script = "scripts/zsh/make_gff.zsh",
        r_script = "scripts/r/trim_gff.R",
        bed = "raw_files/annotations/THEVannotated_genesOnly.bed"
    output:
        "raw_files/annotations/thev_predicted_genes.gff"
    conda:
        "environment.yml"
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
    conda:
        "environment.yml"
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
    conda:
        "environment.yml"
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
    conda:
        "environment.yml"
    shell:
        """
        {input.r_script}
        {input.bash_script}
        """

# #################### BUILD THEV GENOMIC INDEX FOR MAPPING WITH HISAT2 ######
rule build_genome_index:
    input:
        script = "scripts/zsh/build_genome_index.zsh",
        ss = "raw_files/annotations/host_thev.ss",
        exon = "raw_files/annotations/host_thev.exons",
        mix_genome = "raw_files/genome_file/mix_genomes.fa"
    output:
        expand("raw_files/host_virus_genome_index/mix_tran.{n}.ht2", \
        n = range(1,9))
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# #################### MAP READS TO THEV GENOME WITH HISAT2 #############
rule map_sort_to_bam:
    input:
        script = "scripts/zsh/map_sort_to_bam.zsh",
        seqidx = rules.build_genome_index.output,
        fordata = expand("trimmedReads/forwardTrims/LCS9132_I_{tp}hrsS{rep}_Clean_Data1_val_1.fq.gz",
        tp = [72, 24, 4], rep = [1, 2, 3]),
        for12 = expand("trimmedReads/forwardTrims/LCS9132_I_12hrsS{rep}_Clean_Data1_val_1.fq.gz",
        rep = [1, 3] ),
        revdata = expand("trimmedReads/reverseTrims/LCS9132_I_{tp}hrsS{rep}_Clean_Data2_val_2.fq.gz",
        tp = [72, 24, 4], rep = [1, 2, 3]),
        rev12 = expand("trimmedReads/reverseTrims/LCS9132_I_12hrsS{rep}_Clean_Data2_val_2.fq.gz",
        rep = [1, 3])
    output:
        expand("results/hisat2/thev_sorted_{time}hrsS{rep}.bam",
        time = ["72", "24", "4"], rep = ["1", "2", "3"]),
        expand("results/hisat2/thev_sorted_12hrsS{rep}.bam",
        rep = ["1", "3"])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# #################### INDEX ALL SORTED .BAM FILES ##################
rule index_bam_files:
    input:
        script = "scripts/zsh/index.zsh",
        bam = rules.map_sort_to_bam.output
    output:
        expand("results/hisat2/thev_sorted_{time}hrsS{rep}.bam.bai",
        time = ["72", "24", "4"], rep = ["1", "2", "3"]),
        expand("results/hisat2/thev_sorted_12hrsS{rep}.bam.bai",
        rep = ["1", "3"])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# #################### SUBSET SORTED BAM FILES FOR THEV ##################
rule filter_thev:
    input:
        script = "scripts/zsh/thev_subset.zsh",
        bam = rules.map_sort_to_bam.output,
        idx = rules.index_bam_files.output
    output:
        expand("results/hisat2/thev_subset_{time}hrsS{rep}.bam",
        time = ["72", "24", "4"], rep = ["1", "2", "3"]),
        expand("results/hisat2/thev_subset_12hrsS{rep}.bam",
        rep = ["1", "3"])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# #################### INDEX ALL SUBSET BAM FILES ##################
rule index_subset_bam:
    input:
        script = "scripts/zsh/subset_index.zsh",
        bam = rules.filter_thev.output
    output:
        expand("results/hisat2/thev_subset_{time}hrsS{rep}.bam.bai",
        time = ["72", "24", "4"], rep = ["1", "2", "3"]),
        expand("results/hisat2/thev_subset_12hrsS{rep}.bam.bai",
        rep = ["1", "3"])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# #################### CONSTRUCT TRANSCRIPTS WITH STRINGTIE ##############
rule make_transcripts:
    input:
        script = "scripts/zsh/assemble_transcripts.zsh",
        gtf = "raw_files/annotations/thev_predicted_genes.gtf",
        bam = rules.filter_thev.output
    output:
        expand("results/stringtie/thev_{time}hrsS{rep}.gtf",
        time = [24, 72], rep = [1, 2, 3]),
        expand("results/stringtie/thev_12hrsS{rep}.gtf",
        rep = [1, 3]),
        expand("results/stringtie/thev_4hrsS{rep}.gtf",
        rep = [1, 2])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# #################### MERGE ALL GTF FILES #####################
rule merge_gtfs:
    input:
        rules.make_transcripts.output
    output:
        "results/stringtie/all_merged.gtf"
    conda:
        "environment.yml"
    shell:
        "cat {input} > {output}"

# ##################### FILTER FOR REAL TRANSCRIPTS (REMOVE PREDICTED ORFs) ######
rule split_timepoint_transcripts:
    input:
        r_script = "scripts/r/filter_real_transcripts.R",
        all_gtf = rules.merge_gtfs.output
    output:
        "results/stringtie/all_real_transcripts_merged.gtf",
        expand("results/stringtie/transcripts_merged_{tp}hrs.gtf", \
        tp = [4, 12, 24, 72])
    conda:
        "environment.yml"
    shell:
        "{input.r_script}"

# ############## GFFCOMPARE TO GENERATE FINAL TRANSCRIPTOME ##############
rule make_unredundant_transcriptome:
    input:
        gtfs = expand("results/stringtie/transcripts_merged_{tp}hrs.gtf",
        tp = [4, 12, 24, 72]),
        script = "scripts/zsh/gffcompare.zsh",
        orfs = "raw_files/annotations/thev_predicted_genes.gtf"
    output:
        expand("results/gffcompare/gffcomp_alltimes.{metric}",
        metric = ["combined.gtf", "loci", "stats", "tracking"])
    conda:
        "environment.yml"
    shell:
        """
        {input.script}
        mv results/stringtie/gffcomp_alltimes.* results/gffcompare/
        """

# ############## ADD REGION ATTRIBUTE TO GFFCOMPARE GTF FILE ##############
rule mod_final_trxptome:
    input:
        gtf = "results/gffcompare/gffcomp_alltimes.combined.gtf",
        rscript = "scripts/r/add_region_to_trxptome.R"
    output:
        "results/gffcompare/updated_alltimes.combined.gtf"
    conda:
        "environment.yml"
    shell:
        "{input.rscript}"

# ################### MAKE FIGURE 3 (TRANSCRIPTOME MAP) #######################################
rule make_figure3:
    input:
        rscript1 = "scripts/r/plot_thev_timepoint_splicing.R",
        r_pack1 = "scripts/r/geom_brace.R",
        r_pack2 = "scripts/r/seekBrace.R",
        gtfs = expand("results/stringtie/transcripts_merged_{tp}hrs.gtf",
        tp = [4, 12, 24, 72]),
        spliced_gtf = "results/gffcompare/gffcomp_alltimes.combined.gtf",
        orf_gtf = rules.merge_gtfs.output,
        rscript2 = "scripts/r/thev_splicing_fullMap.R"
    output:
        "results/r/figures/figure3.png"
    conda:
        "environment.yml"
    shell:
        "{input.rscript1}"

# ############## ESTIMATE TRANSCRIPT ABUNDANCES##############
rule est_abundances:
    input:
        script = "scripts/zsh/est_trxpt_abund.zsh",
        bams = rules.filter_thev.output,
        ref = rules.mod_final_trxptome.output,
    output:
        expand("results/ballgown/abund_{s}/abund_{s}.gtf",
        s = ["72hrsS1", "72hrsS2", "72hrsS3", "24hrsS1", "24hrsS2", "24hrsS3",
        "12hrsS1", "12hrsS3", "4hrsS1", "4hrsS2"])
    conda:
        "environment.yml"
    shell:
        "{input.script}"


# #################### COUNT ALL SPLICE JUNCTIONS ####################
rule count_junctions:
    input:
        bams = rules.filter_thev.output,
        index = rules.index_subset_bam.output,
        script = "scripts/zsh/splice_site_stats.zsh"
    output:
        expand("results/hisat2/junction_stats_{tp}S{rep}.bed",
        tp = [4, 24, 72], rep = [1, 2, 3]),
        expand("results/hisat2/junction_stats_12S{rep}.bed", rep = [1, 3])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# # #################### BULK MAP READS ##########
rule bulk_map_sort_to_bam:
    input:
        script = "scripts/zsh/bulk_map_sort_to_bam.zsh",
        seqidx = rules.build_genome_index.output,
        fordata = expand("trimmedReads/forwardTrims/LCS9132_I_{tp}hrsS{rep}_Clean_Data1_val_1.fq.gz",
        tp = [72, 24, 4], rep = [1, 2, 3]),
        for12 = expand("trimmedReads/forwardTrims/LCS9132_I_12hrsS{rep}_Clean_Data1_val_1.fq.gz",
        rep = [1, 3]),
        revdata = expand("trimmedReads/reverseTrims/LCS9132_I_{tp}hrsS{rep}_Clean_Data2_val_2.fq.gz",
        tp = [72, 24, 4], rep = [1, 2, 3]),
        rev12 = expand("trimmedReads/reverseTrims/LCS9132_I_12hrsS{rep}_Clean_Data2_val_2.fq.gz",
        rep = [1, 3])
    output:
        expand("results/hisat2/bulk/sortedTHEV_{time}hrsSamples.bam",
        time = [4, 12, 24, 72])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# # #################### BULK INDEX ############
rule bulk_index:
    input:
        script = "scripts/zsh/bulk_index.zsh",
        bam = rules.bulk_map_sort_to_bam.output
    output:
        expand("results/hisat2/bulk/sortedTHEV_{time}hrsSamples.bam.bai",
        time = [4, 12, 24, 72])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# # #################### HOST AND VIRUS BULK COVERAGE #########
rule total_bulk_coverage:
    input:
        script = "scripts/zsh/total_bulk_coverage.zsh",
        bam = rules.bulk_map_sort_to_bam.output,
        idx = rules.bulk_index.output
    output:
        expand("results/hisat2/coverage/host_thev_cov{t}.txt",
        t = [4, 12, 24, 72])
    conda:
        "environment.yml"  
    shell:
        "{input.script}"

# # #################### BULK SUBSET THEV #########
rule filter_bulk_thev:
    input:
        script = "scripts/zsh/subset_bulk_thev.zsh",
        bams = rules.bulk_map_sort_to_bam.output,
        idx = rules.bulk_index.output
    output:
        expand("results/hisat2/bulk/subsetTHEV_{time}hrs.bam",
        time = [4, 12, 24, 72])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# # #################### BULK SUBSET INDEX ############
rule bulk_subset_index:
    input:
        script = "scripts/zsh/subset_bulk_idx.zsh",
        bams = rules.filter_bulk_thev.output
    output:
        expand("results/hisat2/bulk/subsetTHEV_{time}hrs.bam.bai",
        time = [4, 12, 24, 72])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# # #################### BULK DEPTH ############
rule bulk_depth:
    input:
        script = "scripts/zsh/bulk_depth.zsh",
        bam = rules.filter_bulk_thev.output,
        idx = rules.bulk_subset_index.output
    output:
        expand("results/hisat2/coverage/thev_{time}hrsdepth.txt",
        time = [4, 12, 24, 72])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# # #################### BULK COUNT JUNCTIONS ############
rule bulk_count_junctions:
    input:
        bams = rules.filter_bulk_thev.output,
        index = rules.bulk_subset_index.output,
        script = "scripts/zsh/bulk_splice_stats.zsh"
    output:
        expand("results/hisat2/bulk/junction_stats_{tp}hrs.bed",
        tp = [4, 12, 24, 72])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# # #################### BULK ANNOTATE JUNCTIONS ############
rule bulk_annotate_junctions:
    input:
        script = "scripts/zsh/annotate_junctions.zsh",
        bedfiles = rules.bulk_count_junctions.output,
        fasta = "raw_files/genome_file/AY849321.1.fa",
        gtf= rules.mod_final_trxptome.output
    output:
       expand("results/hisat2/bulk/annot_{tp}hrsSS.txt", tp = [4, 12, 24, 72])
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# # #################### BULK COUNTING READS ############
rule count_total_reads:
    input:
        script = "scripts/zsh/count_total_reads.zsh",
        bam = rules.bulk_map_sort_to_bam.output
    output:
        "results/hisat2/coverage/bulk_counts.txt"
    conda:
        "environment.yml"
    shell:
        "{input.script}"

# # #################### PLOT THEV GENOMIC MAP ############
rule make_orf_map:
    input:
        r_script = "scripts/r/thev_orf_map.R",
        bedfile = "raw_files/annotations/THEVannotated_genesOnly.bed",
    output:
        "results/r/figures/figure1.png"
    conda:
        "environment.yml"
    shell:
       "{input.r_script}"

# # #################### PLOT TRANSCRIPT  AND JUNCTION ABUNDANCES OVER TIME ############
rule make_fig4_and_5:
    input:
        trxptome = rules.mod_final_trxptome.output,
        annot_juncs = rules.bulk_annotate_junctions.output,
        abun_files = rules.est_abundances.output,
        bams = rules.filter_thev.output,
        rscript1 = "scripts/r/plot_abundances.R",
        rscript2 = "scripts/r/abundance_analyses.R"
    output:
        expand("results/r/figures/figure{num}.png", num = [4, 5])
    conda:
        "environment.yml"
    shell:
        "{input.rscript1}"

# # #################### MAKE THEV GROWTH CURVE AND MAPPING DEPTH ############
rule make_fig2:
    input:
        r_script = "scripts/r/thev_growthcurve.R",
        exp1 = "raw_files/wetlab_data/thev_growthcurve2021.xls",
        exp2 = "raw_files/wetlab_data/thev_growthcurve04_2023.xls",
        r_script1 = "scripts/r/thev_cov_depth.R",
        r_script2 = "scripts/r/thev_genomic_map.R",
        bedfile = "raw_files/annotations/THEVannotated_genesOnly.bed",
        depth = rules.bulk_depth.output,
        coverage = rules.total_bulk_coverage.output
    output:
        "results/r/figures/figure2.png"
    conda:
        "environment.yml"
    shell:
       "{input.r_script}"

# # #################### PLOT REGION-BY-REGION TRANCRIPTS ############
rule make_fig6_10:
    input:
        r_script = "scripts/r/plot_fig6_10.R",
        r_script1 = "scripts/r/abundance_analyses.R",
        r_script_2 = "scripts/r/reg_by_reg_plots.R",
        trxpts = rules.mod_final_trxptome.output
    output:
        expand("results/r/figures/figure{num}.png",
        num = range(6, 11))
    conda:
        "environment.yml"
    shell:
       "{input.r_script}"

rule make_project_DAG:
    output:
        "project_map.png"
    conda:
        "environment.yml"
    shell:
        """
        snakemake --rulegraph write_manuscript | dot -Tpng > {output}
        """

# # ############### WRITE MANUSCRIPT FOR PUBLICATION ##################
rule write_manuscript:
    input:
        rmd = "manuscript_thev_transcriptome.Rmd",
        ref_style = "asm.csl",
        refs = "transcriptome_refs.bib",
        r1 = "scripts/r/bam_file_analysis.R",
        r2 = "scripts/r/geom_brace.R",
        rbrace_script = "scripts/r/seekBrace.R",
        juncs = rules.bulk_annotate_junctions.output,
        covs = rules.total_bulk_coverage.output,
        tot_reads = rules.count_total_reads.output,
        fig2 = rules.make_fig2.output,
        fig1 = rules.make_orf_map.output,
        fig3 = rules.make_figure3.output,
        fig4_5 = rules.make_fig4_and_5.output,
        abunds = "scripts/r/abundance_analyses.R",
        r3 = "scripts/r/reg_by_reg_plots.R",
        fig6_10 = rules.make_fig6_10.output,
        stylesheet = "style_word_output.docx",
        dag = rules.make_project_DAG.output
    output:
        "manuscript_thev_transcriptome.pdf",
        "manuscript_thev_transcriptome.docx",
        "manuscript_thev_transcriptome.tex"
    conda:
        "environment.yml"
    shell:
        """
        R -e "library(rmarkdown);render('{input.rmd}', output_format = 'all')";
        rm ./manuscript_thev_transcriptome.log
        """

rule write_supplementary:
    input:
        rmd = "supplementary_thev_trxptome.Rmd",
        rscript1 = "scripts/r/bam_file_analysis.R",
        rscript2 = "scripts/r/reg_by_reg_plots.R",
        gels = expand("wet_lab_validation/validation_gels/trxpt_{trx_n}_gel.png",
        trx_n = [1, 2, 3, 5, "6or7", 8, "10or9_j1", 12, 13, "14or11j1", 15,
        16, 17,18, 20, 21, "22or10j2", "23or29", "24or11j2", 25, 26, 27, 28]),
        dag = rules.make_project_DAG.output,
        stylesheet = "style_word_output.docx"
    output:
        "supplementary_thev_trxptome.pdf",
        "supplementary_thev_trxptome.docx"
    conda:
        "environment.yml"
    shell:
        """
        R -e "library(rmarkdown);render('{input.rmd}', output_format = 'all')"
        """

rule study_importance:
    input:
        "importance_of_study.Rmd"
    output:
        "importance_of_study.pdf",
        "importance_of_study.docx"
    conda:
        "environment.yml"
    shell:
        """
        R -e "library(rmarkdown); render('{input}', output_format = 'all')"
        """