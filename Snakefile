import os
from pathlib import Path

configfile: "config.json"

wildcard_constraints:
	var = "({})".format("|".join(config["var"])),
	treat = "({})".format("|".join(config["treat"])),
	day = "({})".format("|".join(config["dai"])),
	N = "({})".format("|".join(config["N"])),
	strand = "({})".format("|".join(config["strand"])),
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
		fastqc = get_fastqc(config)

		quant = expand(
			"results/kallisto/{var}.{treat}.{dai}.{N}.tsv",
			var = config["cultivars"],
			treat = config["treatment"],
			dai = config["dai"],
			N = config["replicates"],
			path = config["paths"],
			format = config["formats"]
		),

		# DESeq2_results_directory = config['path']['DESeq2_results'],
