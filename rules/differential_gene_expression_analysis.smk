

"""
        Differential gene expression analysis
"""


rule DESeq2_init:
    input:
        # combined = rules.combine_gene_quantification.output.combined
        combined = 'data/all.csv'
    output:
        DESeq2_dds_init = config['path']['DESeq2_dds_init'],
        DESeq2_sample_information = config['path']['DESeq2_sample_information']
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_init.R'


rule DESeq2_execute:
    input:
        DESeq2_dds_init = rules.DESeq2_init.output.DESeq2_dds_init
    output:
        DESeq2_dds_execute = config['path']['DESeq2_dds_execute']
    threads:
        6
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_execute.R'


rule DESeq2_get_results:
    input:
        DESeq2_dds_execute = rules.DESeq2_execute.output.DESeq2_dds_execute
    output:
        DESeq2_results_directory = config['path']['DESeq2_results']
    threads:
        6
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_get_results.R'