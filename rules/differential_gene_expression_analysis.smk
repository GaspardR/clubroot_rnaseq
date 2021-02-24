

"""
        Differential gene expression analysis
"""


rule DESeq2_init:
    input:
        #combined = rules.combine_gene_quantification.output.combined
        combined = 'data/all.csv'
    output:
        #DESeq2_dds_init_dir = directory(config['path']['DESeq2_dds_init_dir']),
        sample_information = config['path']['DESeq2_sample_information'],
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
        ),
        dds_init_condition_all_infected = os.path.join(
            config['path']['DESeq2_dds_init_dir'],
            "dds_init_condition_all_infected.dds"
        )
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_init.R'


condition_dai = [
    'all',
    'all_infected',
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

condition_dai = [
    'all',
    '7',
    '14',
    '21'
]

rule DESeq2_get_results:
    input:
        dds_execute = rules.DESeq2_execute.output.dds_execute
    output:
        results = os.path.join(
            config['path']['DESeq2_results'],
            'results_condition_{condition_dai}.csv'
        ),
    threads:
        1
    params:
        condition = '{condition_dai}'
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_get_results.R'

rule DESeq2_get_results_dai:
    input:
        dds_execute = os.path.join(
            config['path']['DESeq2_dds_execute_dir'],
            'dds_execute_condition_all_infected.dds'
        ),
    output:
        results_dai_7vs14 = os.path.join(
            config['path']['DESeq2_results'],
            'results_condition_all_infected_7vs14.csv'
        ),
        results_dai_7vs21 = os.path.join(
            config['path']['DESeq2_results'],
            'results_condition_all_infected_7vs21.csv'
        ),
        results_dai_14vs21 = os.path.join(
            config['path']['DESeq2_results'],
            'results_condition_all_infected_14vs21.csv'
        ),
    threads:
        1
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_get_results_dai.R'

condition_dai = [
    'all',
    '7',
    '14',
    '21',
    'all_infected_7vs14',
    'all_infected_7vs21',
    'all_infected_14vs21'
]

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


rule generate_expression_profile:
    input:
        dds_all = os.path.join(
            config['path']['DESeq2_dds_execute_dir'],
            'dds_execute_condition_all.dds'
        ),
        sample_information = rules.DESeq2_init.output.sample_information,
        results_dai_7vs14 = rules.DESeq2_get_results_dai.output.results_dai_7vs14,
        results_dai_7vs21 = rules.DESeq2_get_results_dai.output.results_dai_7vs21,
        results_dai_14vs21 = rules.DESeq2_get_results_dai.output.results_dai_14vs21,
    output:
        expression_profile = "data/figures/expression_profile.png",
    conda:
        "../envs/DESeq2.yaml"
    script:
        '../modules/R_scripts/generate_expression_profile.R'
