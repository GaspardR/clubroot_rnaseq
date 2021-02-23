

"""
        Differential gene expression analysis
"""


rule DESeq2_init:
    input:
        combined = rules.combine_gene_quantification.output.combined
        #combined = 'data/all.csv'
    output:
        DESeq2_dds_init_dir = directory(config['path']['DESeq2_dds_init_dir']),
        dds_init_condition_all = os.path.join(
            config['path']['DESeq2_dds_init_dir'],
            "dds_init_condition_all.dds"
        ),
        dds_init_condition_7 = os.path.join(
            config['path']['DESeq2_dds_init_dir'],
            "dds_init_condition_7.dds"
        ),
        dds_init_condition_14 = os.path.join(
            config['path']['DESeq2_dds_init_dir'],
            "dds_init_condition_14.dds"
        ),
        dds_init_condition_21 = os.path.join(
            config['path']['DESeq2_dds_init_dir'],
            "dds_init_condition_21.dds"
        )
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_init.R'


condition_dai = [
    'all',
    #'all_infected',
    '7',
    '14',
    '21'
]


rule DESeq2_execute:
    input:
        dds_init = os.path.join(
            config['path']['DESeq2_dds_init_dir'],
            'dds_init_condition_{condition_dai}.dds'
        ),
    output:
        dds_execute = os.path.join(
            config['path']['DESeq2_dds_execute_dir'],
            'dds_execute_condition_{condition_dai}.dds'
        ),
    threads:
        1
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_execute.R'

rule DESeq2_get_results:
    input:
        dds_execute = rules.DESeq2_execute.output.dds_execute
    output:
        results = os.path.join(
            config['path']['DESeq2_results'],
            'results_condition_{condition_dai}.csv'
        ),
    params: 
        condition = '{condition_dai}'
    threads:
        1
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_get_results.R'

rule generate_volcanoplot:
    input:
        results = rules.DESeq2_get_results.output.results
    output:
        volcanoplot = os.path.join(
            config['path']['volcanoplots'],
            'volcanoplot_condition_{condition_dai}.png'
        ),

    conda:
        "../envs/DESeq2.yaml"
    script:
        '../modules/R_scripts/generate_volcanoplot.R'