from pathlib import Path


rule download_sra:
	output:
		fq = Path(config["path"]["raw"], "{var}.{treat}.{dai}.{N}.1")
	run:
		key = '.'.join(
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

"""
	Download the genomes of each specie
"""
rule download_genome_brassica_napus:
    output:
        genome_brassica_napus = config["path"]["genome_brassica_napus"]
    shell:
        "wget --quiet -O {output.genome_brassica_napus}.gz {config[download][genome_brassica_napus]}"
        " && gunzip {output.genome_brassica_napus}.gz"

rule download_genome_plasmodiophora_brassicae:
    output:
        genome_plasmodiophora_brassicae = config["path"]["genome_plasmodiophora_brassicae"]
    shell:
        "wget --quiet -O {output.genome_plasmodiophora_brassicae}.gz {config[download][genome_plasmodiophora_brassicae]}"
        " && gunzip {output.genome_plasmodiophora_brassicae}.gz"

"""
	Download the annotations of each specie
"""
rule download_annotation_brassica_napus:
	output:
		gtf_brassica_napus = config["path"]["gtf_brassica_napus"]
	shell:
		"wget --quiet -O {output.gtf_brassica_napus}.gz {config[download][gtf_brassica_napus]}"
        " && gunzip {output.gtf_brassica_napus}.gz"

rule download_annotation_plasmodiophora_brassicae:
	output:
		gtf_plasmodiophora_brassicae = config["path"]["gtf_plasmodiophora_brassicae"],
		gtf_plasmodiophora_brassicae_noheader = config['path']['gtf_plasmodiophora_brassicae_noheader']
	shell:
		"wget --quiet -O {output.gtf_plasmodiophora_brassicae}.gz {config[download][gtf_plasmodiophora_brassicae]}"
        " && gunzip {output.gtf_plasmodiophora_brassicae}.gz"
		" && grep -v '#!' {output.gtf_plasmodiophora_brassicae} > {output.gtf_plasmodiophora_brassicae_noheader}"


"""
	Merge the annotation and the genome together
"""

rule merge_annotation:
	input:
		annotation_brassica_napus = rules.download_annotation_brassica_napus.output.gtf_brassica_napus,
		annotation_plasmodiophora_brassicae = rules.download_annotation_plasmodiophora_brassicae.output.gtf_plasmodiophora_brassicae_noheader
	output:
		merged_annotation = config['path']['merged_annotation']
	shell:
		"cat {input.annotation_plasmodiophora_brassicae} {input.annotation_brassica_napus} > {output.merged_annotation}"


rule merge_genome:
	input:
		genome_plasmodiophora_brassicae = rules.download_genome_plasmodiophora_brassicae.output.genome_plasmodiophora_brassicae,
		genome_brassica_napus = rules.download_genome_brassica_napus.output.genome_brassica_napus
	output:
		merged_genome = config['path']['merged_genome']
	shell:
		"cat {input.genome_brassica_napus} {input.genome_plasmodiophora_brassicae} > {output.merged_genome}"

