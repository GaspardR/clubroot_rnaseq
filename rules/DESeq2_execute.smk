

"""
        Differential gene expression analysis using the DESeq2 R package. Generate the "dds" object.
"""


rule DESeq2_execute:
    input:
        DESeq2_dds_directory = rules.DESeq2_init.output.DESeq2_dds_directory
    output:
        DESeq2_results_directory = "data/DESeq2/results"
    params:
        cores = 3
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_execute.R'