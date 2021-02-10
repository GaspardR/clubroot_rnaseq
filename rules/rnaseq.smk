import os
from pathlib import Path


rule trimmomatic:
    input:
        r = rules.fasterq_dump.output
    output:
        r = Path(
            config["path"]["trimmed"],
            "{var}.{treat}.{dai}.{N}.trim.fastq"
        )
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "logs/trimmomatic/{var}.{treat}.{dai}.{N}.log"
    params:
        trimmer = [
            "ILLUMINACLIP:{}:2:30:10".format(config['path']['adapters']),
            "LEADING:30",
            "TRAILING:30",
            "SLIDINGWINDOW:5:30",
            "MINLEN:15"
        ]
    threads:
        8
    shell:
        "trimmomatic SE"
        " -threads {threads}"
        " -phred33"
        " {input.r}"
        " {output.r}"
        " {params.trimmer}"
        " &> {log}"


rule fastqc:
    input:
        fq = rules.fasterq_dump.output,
        trim = rules.trimmomatic.output.r
    output:
        raw_html = Path(
            config["path"]["qc"],
            "{var}.{treat}.{dai}.{N}_1_fastqc.html"
        ),
        raw_zip = Path(
            config["path"]["qc"],
            "{var}.{treat}.{dai}.{N}_1_fastqc.zip"
        ),
        trim_html = Path(
            config["path"]["qc"],
            "{var}.{treat}.{dai}.{N}.trim_fastqc.html"
        ),
        trim_zip = Path(
            config["path"]["qc"],
            "{var}.{treat}.{dai}.{N}.trim_fastqc.zip"
        ),
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/{var}.{treat}.{dai}.{N}.log"
    params:
        outdir = config["path"]["qc"]
    threads:
        8
    shell:
        "fastqc"
        " --outdir {params.outdir}"
        " --format fastq"
        " --threads {threads}"
        " {input.fq} {input.trim}"
        " &> log"


rule create_transcriptome:
    """ Uses gffread to generate a transcriptome """
    input:
        genome = rules.merge_genome.output.merged_genome,
        gtf = rules.merge_annotation.output.merged_annotation
    output:
        transcriptome = config["path"]["transcriptome"]
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output.transcriptome}"


rule kallisto_index:
    """ Generates the transcriptome index for Kallisto """
    input:
        transcriptome = rules.create_transcriptome.output.transcriptome
    output:
        kallisto_idx = config['path']['kallisto_idx']
    params:
        kmer = "31"
    log:
        "logs/kallisto/index.log"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index "
        "--index={output.kallisto_idx} "
        "--kmer-size={params.kmer} "
        "{input.transcriptome} "
        "&> {log}"


rule kallisto_quant:
    """ Generates counts using Kallisto pseudo-alignment """
    input:
        kallisto_idx = rules.kallisto_index.output.kallisto_idx,
        fq = rules.trimmomatic.output.r
    output:
        quant = Path(
            config["path"]["kallisto_quant"],
            "{var}.{treat}.{dai}.{N}/abundance.tsv"
        ),
        h5 = Path(
            config["path"]["kallisto_quant"],
            "{var}.{treat}.{dai}.{N}/abundance.h5"
        )
    params:
        bootstrap = "50",
        outdir = Path(
            config["path"]["kallisto_quant"],
            "{var}.{treat}.{dai}.{N}"
        )
    log:
        "logs/kallisto/{var}.{treat}.{dai}.{N}.log"
    threads:
        config["params"]["kallisto_quant"]["cores"]
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.kallisto_idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "--single -l 200 -s 20 "
        "{input.fq} "
        "&> {log}"


rule combine_gene_quantification:
    input:
        counts = expand(
            os.path.join(
                config["path"]["kallisto_quant"],
                "{var}.{treat}.{dai}.{N}/abundance.tsv"
            ),
            var = config["var"],
            treat = config["treat"],
            dai = config["dai"],
            N = config["N"]
        ),
        map = rules.transcriptID2geneName.output.map
    output:
        combined = config["path"]["combined"]
    conda:
        "../envs/pypackages.yaml"
    script:
        "../modules/python_scripts/combine_gene_quantification.py"
