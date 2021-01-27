import os
from pathlib import Path

configfile: "config.json"

wildcard_constraints:
	var = "({})".format("|".join(config["cultivars"])),
	treat = "({})".format("|".join(config["treatment"])),
	day = "({})".format("|".join(config["dai"])),
	N = "({})".format("|".join(config["replicates"])),
	strand = "({})".format("|".join(config["strands"])),
	path = "({})".format("|".join(config["paths"])),
	format = "({})".format("|".join(config["formats"]))


include: "rules/downloads.smk"
include: "rules/rnaseq.smk"
include: "rules/differential_gene_expression_analysis.smk"


__author__ = ["Gaspard Reulet", "Hoang-Dong Nguyen"]


rule all:
	input:
		fastqc = expand(
			os.path.join(
				config["path"]["qc"],
	            "{var}.{treat}.{dai}.{N}.{path}_fastqc.{format}"
			),
			var = config["cultivars"],
			treat = config["treatment"],
			dai = config["dai"],
			N = config["replicates"],
			path = config["paths"],
			format = config["formats"]
		),

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


