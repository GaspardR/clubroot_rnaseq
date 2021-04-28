

"""
        Differential gene expression analysis
"""


rule DESeq2_init:
    input:
        combined = 'data/all.csv'
    output:
        DESeq2_dds_init_dir = directory(config['path']['DESeq2_dds_init_dir']),
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
        DESeq2_size_factor_normalization = os.path.join(
            config['path']['DESeq2_size_factor_normalization'],
            'size_factor_normalization_count_{condition_dai}.csv'
        ),
        DESeq2_vst_normalization = os.path.join(
            config['path']['DESeq2_vst_normalization'],
            'vst_normalization_count_{condition_dai}.csv'
        ),

    params:
        condition = '{condition_dai}'
    threads:
        1
    conda:
        '../envs/DESeq2.yaml'
    script:
        '../modules/R_scripts/DESeq2_execute.R'

#condition_dai = [
#    'all',
#    '7',
#    '14',
#    '21'
#]

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

#condition_dai2 = [
#    'all',
#    '7',
#    '14',
#    '21',
#    'all_infected_7vs14',
#    'all_infected_7vs21',
#    'all_infected_14vs21'
#]

#rule generate_volcanoplot:
#    input:
#        results = rules.DESeq2_get_results.output.results
#    output:
#        volcanoplot = os.path.join(
#            config['path']['volcano_plots'],
#            'volcanoplot_condition_{condition_dai2}.png'
#        ),
#    conda:
#        "../envs/DESeq2.yaml"
#    script:
#        '../modules/R_scripts/generate_volcanoplot.R'


rule generate_expression_profile:
    input:
        #sizefactor_normalized_count = 'data/DESeq2/DESeq2_sizefactor_normalized_count.csv',
        sample_information = rules.DESeq2_init.output.sample_information,
        results_dai_7vs14 = rules.DESeq2_get_results_dai.output.results_dai_7vs14,
        results_dai_7vs21 = rules.DESeq2_get_results_dai.output.results_dai_7vs21,
        results_dai_14vs21 = rules.DESeq2_get_results_dai.output.results_dai_14vs21,
        sizefactor_normalized_count = os.path.join(
            config['path']['DESeq2_size_factor_normalization'],
            'size_factor_normalization_count_all.csv'
        ),
    output:
        expression_profile = directory(
            config['path']['expressionprofil_plots']
        ),
    conda:
        "../envs/DESeq2.yaml"
    script:
        '../modules/R_scripts/generate_expression_profile.R'

rule generate_heatmap:
    input:
        generate_expression_dir = rules.generate_expression_profile.output.expression_profile,
        #vst_normalized_count = 'data/DESeq2/DESeq2_vst_normalized_count.csv',
        vst_normalized_count = os.path.join(
            config['path']['DESeq2_vst_normalization'],
            'vst_normalization_count_all.csv'
        ),
    output:
        heatmap_dir = directory(
            config['path']['heatmap_dir']
        ),
    conda:
        "../envs/heatmap.yaml"
    script:
        '../modules/R_scripts/generate_heatmap.R'

