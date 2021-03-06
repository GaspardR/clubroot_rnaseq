import os
from pathlib import Path

configfile: "config.json"

wildcard_constraints:
	var = "({})".format("|".join(config["var"])),
	treat = "({})".format("|".join(config["treat"])),
	dai = "({})".format("|".join(config["dai"])),
	N = "({})".format("|".join(config["N"])),
	format = "({})".format("|".join(config["format"]))


include: "rules/downloads.smk"
include: "rules/rnaseq.smk"
include: "rules/differential_gene_expression_analysis.smk"


__author__ = ["Gaspard Reulet", "Hoang-Dong Nguyen"]


def get_fastqc(config):
	files = list()
	path = config["path"]["qc"]
	raw = "/{var}.{treat}.{dai}.{N}_1_fastqc.{format}"
	trim = "/{var}.{treat}.{dai}.{N}.trim_fastqc.{format}"
	files.extend(expand(path + raw, **config))
	files.extend(expand(path + trim, **config))
	return files


rule all:
	input:
		volcanoplot = expand(
				os.path.join(
            	config['path']['volcanoplots'],
            	'volcanoplot_condition_{condition_dai}.png'
        	),
			condition_dai = condition_dai
		),
		fastqc = get_fastqc(config),
		quant = expand(
			os.path.join(
				config["path"]["kallisto_quant"],
				"{var}.{treat}.{dai}.{N}/abundance.tsv"
			),
			var = config["var"],
			treat = config["treat"],
			dai = config["dai"],
			N = config["N"],
		),
		combined = config["path"]["combined"],
		bedgraph = expand(
			os.path.join(
				config["path"]["bedgraph"],
				"{var}.{treat}.{dai}.{N}.bedgraph"
			),
			var = config["var"],
			treat = config["treat"],
			dai = config["dai"],
			N = config["N"],
		),
