from pathlib import Path


rule download_sra:
	output:
		fq = Path(config["path"]["raw"], "{var}.{treat}.{day}.{N}")
	run:
		key = '.'.join(
			[wildcards.var, wildcards.treat, wildcards.day, wildcards.N]
		)
		link = config["files"][key]
		wget_file = Path(output.fq.parent, f"{key}.1")
		shell(f"wget -0 {wget_file} {link}")


rule fasterq_dump:
	input:
		rules.download_sra.output
	output:
		Path(config["path"]["raw"], "{var}.{treat}.{day}.{N}_1.fastq")
	params:
		dir = config["path"]["raw"],
		cores = config["params"]["fasterq_dump"]["cores"]
	conda:
		"../envs/sra-tools.yaml"
	shell:
		"fasterq_dump --threads {params.cores} --dir {params.dir} {input}"
