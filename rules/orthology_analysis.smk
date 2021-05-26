
"""
        Orthology analysis on the gene of interest
"""

rule create_bed_file:
    input:
        heatmap_dir = rules.generate_heatmap.output.heatmap_dir,
        gtf_plasmodiophora_brassicae = config["path"]["gtf_plasmodiophora_brassicae"]
    output: 
        p_brassicae_bed_file = config["path"]["p_brassicae_bed_file"]
    conda:
        "../envs/heatmap.yaml"
    script:
        "../modules/R_scripts/create_bed_file.R"


rule extract_sequence_from_bed:
    input:
        p_brassicae_bed_file = rules.create_bed_file.output.p_brassicae_bed_file,
        genome_plasmodiophora_brassicae = config["path"]["genome_plasmodiophora_brassicae"]
    output:
        p_brassicae_sequences = config["path"]["p_brassicae_sequences"]
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools getfasta "
        "-fi {input.genome_plasmodiophora_brassicae} "
        "-bed {input.p_brassicae_bed_file} "
        "-name "
        "-tab "
        "-fo {output.p_brassicae_sequences}"