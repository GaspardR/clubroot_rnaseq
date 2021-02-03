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
		DESeq2_results_directory = config['path']['DESeq2_results'],
		bedgraph = expand(
			os.path.join(
	            config["path"]["bedgraph"],
	            "{var}.{treat}.{dai}.{N}.{strand}.bedgraph"
	        ),
			var = config["var"],
			treat = config["treat"],
			dai = config["dai"],
			N = config["N"],
			strand = config["strand"]
		),
		# fastqc = expand(
		# 	os.path.join(
		# 		config["path"]["qc"],
	    #         "{var}.{treat}.{dai}.{N}.{path}_fastqc.{format}"
		# 	),
		# 	var = config["cultivars"],
		# 	treat = config["treatment"],
		# 	dai = config["dai"],
		# 	N = config["replicates"],
		# 	path = config["paths"],
		# 	format = config["formats"]
		# )
		fastqc = get_fastqc(config)
		# raw_html = Path(
        #     config["path"]["qc"],
        #     "{var}.{treat}.{dai}.{N}_1_fastqc.html"
        # ),
        # raw_zip = Path(
        #     config["path"]["qc"],
        #     "{var}.{treat}.{dai}.{N}_1_fastqc.zip"
        # ),
        # trim_html = Path(
        #     config["path"]["qc"],
        #     "{var}.{treat}.{dai}.{N}.trim_fastqc.html"
        # ),
        # trim_zip = Path(
        #     config["path"]["qc"],
        #     "{var}.{treat}.{dai}.{N}.trim_fastqc.zip"
        # ),
