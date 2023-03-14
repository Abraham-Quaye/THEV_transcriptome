rule make_gff:
    input:
        script = "scripts/zsh/make_gff.zsh",
        bedfile = "raw_files/annotations/THEVannotated_genesOnly.bed"
    output:
        "raw_files/annotations/thev_predicted_genes.gff" 
    shell:
        """
        {input.script}
        """

rule make_gtf:
    input:
        script = "scripts/zsh/make_gtf.zsh",
        gff = "raw_files/annotations/thev_predicted_genes.gff"
    output:
        "raw_files/annotations/thev_predicted_genes.gtf" 
    shell:
        """
        {input.script}
        """

rule extract_ss:
    input:
        script = "scripts/zsh/extract_ss.zsh",
        gtf = "raw_files/annotations/thev_predicted_genes.gtf"
    output:
        "raw_files/annotations/thev_predicted_genes.ss"
    shell:
        """
        {input.script}
        """

rule extract_exons:
    input:
        script = "scripts/zsh/extract_exons.zsh",
        gtf = "raw_files/annotations/thev_predicted_genes.gtf"
    output:
        "raw_files/annotations/thev_predicted_genes.exons"
    shell:
        """
        {input.script}
        """

rule build_genome_index:
    input:
        script = "scripts/zsh/build_genome_index.zsh",
        ss = "raw_files/annotations/thev_predicted_genes.ss",
        exons = "raw_files/annotations/thev_predicted_genes.exons",
        genome = "raw_files/genome_file/THEV.fa"
    output:
        "raw_files/thevgenome_index/thev_tran"
    shell:
        """ 
        {input.script}
        """
# f = []
# for file in (trimmedReads/forwardTrims/*val_1.fq.gz):
#     f.append(file)

# r = []
# for file in (trimmedReads/reverseTrims/*val_2.fq.gz):
#     r.append(file)

rule mapping:
    input:
        script = "scripts/zsh/mapping.zsh"
    output:
        "results/hisat2/*.sam"
    shell:
        """
        {input.script}
        """