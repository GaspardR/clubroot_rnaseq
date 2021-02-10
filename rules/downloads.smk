from pathlib import Path


rule download_sra:
	output:
		fq = Path(config["path"]["raw"], "{var}.{treat}.{dai}.{N}.1")
	run:
		key = ".".join(
			[wildcards.var, wildcards.treat, wildcards.dai, wildcards.N]
		)
		link = config["download"]["files"][key]
		wget_file = Path(Path(output.fq).parent, f"{key}.1")
		shell(f"wget -O {wget_file} {link}")


rule fasterq_dump:
	input:
		rules.download_sra.output
	output:
		Path(config["path"]["raw"], "{var}.{treat}.{dai}.{N}_1.fastq")
	params:
		cores = config["params"]["fasterq_dump"]["cores"]
	conda:
		"../envs/sra-tools.yaml"
	shell:
		"fasterq-dump --threads {params.cores} --outfile {output} {input}"


# Download the genomes of each species
rule download_genome_brassica_napus:
	output:
		genome = config["path"]["genome_brassica_napus"]
	params:
		link = config["download"]["genome_brassica_napus"]
	shell:
		"wget --quiet -O {output.genome}.gz {params.link}"
		" && gunzip {output.genome}.gz"

rule download_genome_plasmodiophora_brassicae:
	output:
		genome = config["path"]["genome_plasmodiophora_brassicae"]
	params:
		link = config["download"]["genome_plasmodiophora_brassicae"]
	shell:
		"wget --quiet -O {output.genome}.gz {params.link}"
		" && gunzip {output.genome}.gz"


# Download the annotations of each species
rule download_annotation_brassica_napus:
	output:
		gtf = config["path"]["gtf_brassica_napus"]
	params:
		link = config["download"]["gtf_brassica_napus"]
	shell:
		"wget --quiet -O {output.gtf}.gz {params.link}"
		" && gunzip {output.gtf}.gz"

rule download_annotation_plasmodiophora_brassicae:
	output:
		gtf = config["path"]["gtf_plasmodiophora_brassicae"],
	params:
		link = config["download"]["gtf_plasmodiophora_brassicae"]
	shell:
		"wget --quiet -O {output.gtf}.gz {params.link}"
		" && zgrep -v '#!' {output.gtf}.gz > {output.gtf}"


# Merge annotations and genomes
rule merge_annotation:
	input:
		brassica_napus = rules.download_annotation_brassica_napus.output.gtf,
		plasmodiophora_brassicae = rules.download_annotation_plasmodiophora_brassicae.output.gtf
	output:
		merged_annotation = config["path"]["merged_annotation"]
	shell:
		"cat {input.plasmodiophora_brassicae} {input.brassica_napus}"
		" > {output.merged_annotation}"

rule merge_genome:
	input:
		plasmodiophora_brassicae = rules.download_genome_plasmodiophora_brassicae.output.genome,
		brassica_napus = rules.download_genome_brassica_napus.output.genome
	output:
		merged_genome = config["path"]["merged_genome"]
	shell:
		"cat {input.brassica_napus} {input.plasmodiophora_brassicae}"
		" > {output.merged_genome}"


rule transcriptID2geneName:
	input:
		gtf = config["path"]["merged_annotation"]
	output:
		map = config["path"]["map"]
	conda:
		"../envs/pypackages.yaml"
	script:
		"../modules/python_scripts/transcriptID2geneName.py"
