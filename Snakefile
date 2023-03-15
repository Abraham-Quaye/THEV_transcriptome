import glob
import os


rule make_gff:
    input:
       "scripts/zsh/make_gff.zsh"
    output:
        "raw_files/annotations/thev_predicted_genes.gff"
    params:
        "raw_files/annotations/THEVannotated_genesOnly.bed"
    shell:
        "{input}"

rule make_gtf:
    input:
        script = "scripts/zsh/make_gtf.zsh",
        gff = "raw_files/annotations/thev_predicted_genes.gff"
    output:
        "raw_files/annotations/thev_predicted_genes.gtf"
    shell:
        "{input.script}"

rule extract_splice_site:
    input:
        script = "scripts/zsh/extract_ss.zsh",
        gtf = "raw_files/annotations/thev_predicted_genes.gtf"
    output:
        "raw_files/annotations/thev_predicted_genes.ss"
    shell:
        "{input.script}"

rule extract_exons:
    input:
        script = "scripts/zsh/extract_exons.zsh",
        gtf = "raw_files/annotations/thev_predicted_genes.gtf"
    output:
        "raw_files/annotations/thev_predicted_genes.exons"
    shell:
        "{input.script}"

rule build_index:
    input:
        script = "scripts/zsh/build_genome_index.zsh",
        ss = "raw_files/annotations/thev_predicted_genes.ss",
        exon = "raw_files/annotations/thev_predicted_genes.exons"
    params:
        "raw_files/genome_file/THEV.fa"
    output:
        idx = "raw_files/thevgenome_index/*.ht2"
    shell:
        "{input.script}"

