import os
from pathlib import Path

configfile: "config.json"

wildcard_constraints:
	var = "({})".format("|".join(config["cultivars"])),
	treat = "({})".format("|".join(config["treatment"])),
	day = "({})".format("|".join(config["dai"])),
	N = "({})".format("|".join(config["replicates"]))


include: "rules/downloads.smk"

__author__ = ["Gaspard Reulet", "Hoang-Dong Nguyen"]


rule all:
	input:
		fq = expand(
			os.path.join(config["path"]["raw"], "{var}.{treat}.{dai}.{N}_1.fastq"),
			var = config["cultivars"],
			treat = config["treatment"],
			dai = config["dai"],
			N = config["replicates"]
		)
