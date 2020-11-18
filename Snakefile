from pathlib import Path

configfile: "config.json"

include: "rules/download.smk"

__author__ = ["Gaspard Reulet", "Hoang-Dong Nguyen"]


rule all:
	Path(config["path"]["raw"], "{var}.{treat}.{day}.{N}.1.fq.gz")
