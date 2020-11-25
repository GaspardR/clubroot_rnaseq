

"""
        Differential gene expression analysis using the DESeq2 R package. Generate the "dds" object.
"""


rule DESeq2_init:
    # input:
    #     rawcount_matrix = rules.RULES_NAMES.output.RAWCOUNT_OUTPUTNAMES
    output:
        DESeq2_dds_directory = "data/DESeq2/dds"
    params:
        cores = 3
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_init.R'