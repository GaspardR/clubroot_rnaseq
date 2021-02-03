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


# rule STAR_index:
#     input:
#         genome = rules.download_genome.output,
#         gtf = rules.download_annotation.output
#     output:
#         directory(config['path']['STAR_ref'])
#     conda:
#         "../envs/STAR.yaml"
#     threads:
#         16
#     log:
#         "logs/STAR_index.log"
#     shell:
#         "mkdir {output}"
#         " && STAR"
#         " --runMode genomeGenerate"
#         " --runThreadN {threads}"
#         " --genomeDir {output}"
#         " --genomeFastaFiles {input.genome}"
#         " --sjdbGTFfile {input.gtf}"
#         " --sjdbOverhang 74"


# rule STAR_align:
#     input:
#         r = rules.trimmomatic.output.r,
#         index = Path(config['path']['STAR_ref'])
#     output:
#         bam = Path(
#             config["path"]["STAR_align"],
#             "{var}.{treat}.{dai}.{N}.Aligned.sortedByCoord.out.bam"
#         )
#     conda:
#         "../envs/STAR.yaml"
#     params:
#         outFileNamePrefix = os.path.join(
#             config["path"]["STAR_align"],
#             "{var}.{treat}.{dai}.{N}."
#         )
#     threads:
#         16
#     log:
#         "logs/STAR_align.{var}.{treat}.{dai}.{N}.log"
#     shell:
#         "STAR"
#         " --runMode alignReads"
#         " --genomeDir {input.index}"
#         " --readFilesIn {input.r}"
#         " --outFileNamePrefix {params.outFileNamePrefix}"
#         " --runThreadN {threads}"
#         " --readFilesCommand zcat"
#         " --outReadsUnmapped Fastx"
#         " --outStd Log"
#         " --outSAMtype BAM SortedByCoordinate"
#         " --outSAMunmapped None"
#         " --outFilterType BySJout"
#         " --outFilterMismatchNmax 5"
#         " --alignIntronMax 2500"
#         " --alignMatesGapMax 14750"
#         " &> {log}"


# rule coco_ca:
#     input:
#         coco = rules.download_coco.output,
#         gtf = rules.download_annotation.output
#     output:
#         ref = config["path"]["coco_ca"]
#     conda:
#         '../envs/coco.yaml'
#     shell:
#         'python {input.coco} ca'
#         ' -b snoRNA,tRNA,snRNA,ncRNA'
#         ' -o {output.ref}'
#         ' {input.gtf}'


# rule coco_cc:
#     input:
#         gtf = rules.coco_ca.output.ref,
#         bam = rules.STAR_align.output.bam
#     output:
#         counts = Path(
#             config["path"]["coco_cc"],
#             '{var}.{treat}.{dai}.{N}.tsv'
#         )
#     params:
#         coco_path = 'other_git_repos/coco/bin'
#     threads:
#         8
#     conda:
#         '../envs/coco.yaml'
#     shell:
#         'python {params.coco_path}/coco.py cc'
#         ' --countType both'
#         ' --thread {threads}'
#         ' --strand 1'
#         ' {input.gtf}'
#         ' {input.bam}'
#         ' {output.counts}'


rule create_transcriptome:
    """ Uses gffread to generate a transcriptome """
    input:
        genome = rules.download_genome.output,
        gtf = rules.download_annotation.output
    output:
        transcriptome = config['path']['transcriptome']
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
        quant = "results/kallisto/{var}.{treat}.{dai}.{N}.tsv",
        h5 = "results/kallisto/{var}.{treat}.{dai}.{N}.h5",
    params:
        bootstrap = "50",
        outdir = "results/kallisto"
    log:
        "logs/kallisto/{var}.{treat}.{dai}.{N}.log"
    threads:
        32
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







# rule combine_gene_quantification:
#     input:
#         counts = expand(
#             os.path.join(
#                 config['path']['coco_cc'],
#                 '{var}.{treat}.{dai}.{N}.tsv'
#             ),
#             var = config['cultivars'],
#             treat = config['treatment'],
#             dai = config['dai'],
#             N = config['replicates']
#         )
#     output:
#         combined = config['path']['combined']
#     conda:
#         '../envs/pypackages.yaml'
#     script:
#         '../modules/python_scripts/combine_gene_quantification.py'
