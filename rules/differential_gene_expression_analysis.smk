

"""
        Differential gene expression analysis
"""


rule DESeq2_init:
    input:
        # combined = rules.combine_gene_quantification.output.combined
        combined = 'data/test.csv'
    output:
        DESeq2_dds_directory = config['path']['DESeq2_dds']
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_init.R'


rule DESeq2_execute:
    input:
        DESeq2_dds_directory = rules.DESeq2_init.output.DESeq2_dds_directory
    output:
        DESeq2_results_directory = config['path']['DESeq2_results']
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_execute.R'
