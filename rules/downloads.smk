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


rule download_genome:
    """ Download genome """
    output:
        config["path"]["genome"]
    shell:
        "wget --quiet -O {output}.gz {config[download][genome]}"
        " && gunzip {output}.gz"


rule download_annotation:
	output:
		config["path"]["gtf"]
	shell:
		"wget --quiet -O {output}.gz {config[download][gtf]}"
        " && gunzip {output}.gz"


rule download_coco:
    output:
        config["path"]["coco"]
    params:
        dir = Path(config["path"]["coco"]).parent.parent.parent
    conda:
        "../envs/git.yaml"
    shell:
        "mkdir -p {params.dir}"
        " && cd {params.dir}"
        " && rm -rf coco/"
        " && git clone {config[download][coco]}"
