
"""
        Orthology analysis on the gene of interest
"""

rule create_bed_file:
    input:
        heatmap_dir = rules.generate_heatmap.output.heatmap_dir,
        gtf_plasmodiophora_brassicae = config["path"]["gtf_plasmodiophora_brassicae"]
    output: 
        bed_file = config["path"]["bedfile"]
    conda:
        "../envs/heatmap.yaml"
    script:
        "../modules/R_scripts/create_bed_file.R"

rule extract_sequence_from_bed:
    input:
        bed_file = rules.create_bed_file.output.bed_file,
        genome_plasmodiophora_brassicae = config["path"]["genome_plasmodiophora_brassicae"]
    output:
        clubroot_gene_sequences = "data/clubroot_gene_sequences.tsv"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools getfasta "
        "-fi {input.genome_plasmodiophora_brassicae} "
        "-bed {input.bed_file} "
        "-name "
        "-tab "
        "-fo {output.clubroot_gene_sequences}"