
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


rule extract_go_id_from_genes:
    input:
        #heatmap_dir = rules.generate_heatmap.output.heatmap_dir,
        goa_plasmodiophora_brassicae = rules.download_goa_plasmodiophora_brassicae.output.goa_plasmodiophora_brassicae
    output:
        go_id = config["path"]["go_id"]
    conda:
        "../envs/orthology.yaml"
    script:
        "../modules/R_scripts/extract_go_id_from_genes.R"


rule go_enrichment:
    input:
        heatmap_dir = rules.generate_heatmap.output.heatmap_dir,
        go_id = rules.extract_go_id_from_genes.output.go_id
    output:
        enrichment_results = config["path"]["enrichment_results"]
    conda:
        "../envs/orthology.yaml"
    script:
        "../modules/R_scripts/go_enrichment.R"
