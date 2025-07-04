# Snakefile.qc: Quality control + trimming

import pandas as pd
import os

# CONFIG
# config file for directories

configfile: "config/config.yaml"

RAW_DIR = config["raw_dir"]
TRIMMED_DIR = config["trimmed_dir"]
FASTQC_DIR = config["fastqc_dir"]
TRIMMED_FASTQC_DIR = config["trimmed_fastqc_dir"]
MULTIQC_DIR = config["multiqc_dir"]


# INPUT: Find all raw FastQ files
# sample file for all samples

samples = pd.read_csv(config["samples"], sep="\t")
SAMPLES = samples["sample"].tolist()

# RULE: Define workflow

rule all:
    input:
        expand(os.path.join(FASTQC_DIR, "{sample}_R1_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(FASTQC_DIR, "{sample}_R2_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(FASTQC_DIR, "{sample}_R1_fastqc.zip"), sample=SAMPLES),
        expand(os.path.join(FASTQC_DIR, "{sample}_R2_fastqc.zip"), sample=SAMPLES),
        expand(os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(TRIMMED_FASTQC_DIR, "{sample}_R1.trimmed_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(TRIMMED_FASTQC_DIR, "{sample}_R2.trimmed_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(TRIMMED_FASTQC_DIR, "{sample}_R1.trimmed_fastqc.zip"), sample=SAMPLES),
        expand(os.path.join(TRIMMED_FASTQC_DIR, "{sample}_R2.trimmed_fastqc.zip"), sample=SAMPLES),
        os.path.join(MULTIQC_DIR, "multiqc_report.html")


# RULE: FastQC on raw reads

rule fastqc_raw:
    input:
        fq = lambda wildcards: samples.loc[samples["sample"] == wildcards.sample, f"{wildcards.read}"].values[0]
    output:
        html = FASTQC_DIR + "/{sample}_{read}_fastqc.html",
        zip = FASTQC_DIR + "/{sample}_{read}_fastqc.zip"
    shell:
        """
        fastqc {input.fq} --outdir {FASTQC_DIR}
        """

# RULE: Trim adapters with fastp

rule trim_reads:
    input:
        R1 = lambda wildcards: samples.loc[samples["sample"] == wildcards.sample, "R1"].values[0],
        R2 = lambda wildcards: samples.loc[samples["sample"] == wildcards.sample, "R2"].values[0]
    output:
        R1_trimmed = TRIMMED_DIR + "/{sample}_R1.trimmed.fastq.gz",
        R2_trimmed = TRIMMED_DIR + "/{sample}_R2.trimmed.fastq.gz",
        html = TRIMMED_DIR + "/{sample}_fastp.html",
        json = TRIMMED_DIR + "/{sample}_fastp.json"
    threads: 4
    shell:
        """
        fastp \
          -i {input.R1} -I {input.R2} \
          -o {output.R1_trimmed} -O {output.R2_trimmed} \
          -h {output.html} -j {output.json} \
          -w {threads}
        """


# RULE: FastQC on trimmed reads

rule fastqc_trimmed:
    input:
        fq = TRIMMED_DIR + "/{sample}_{read}.trimmed.fastq.gz"
    output:
        html = TRIMMED_FASTQC_DIR + "/{sample}_{read}.trimmed_fastqc.html",
        zip = TRIMMED_FASTQC_DIR + "/{sample}_{read}.trimmed_fastqc.zip"
    shell:
        """
        fastqc {input.fq} --outdir {TRIMMED_FASTQC_DIR}
        """

# RULE: MultiQC summary

rule multiqc:
    input:
        raw = expand(FASTQC_DIR + "/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
        trimmed = expand(TRIMMED_FASTQC_DIR + "/{sample}_{read}.trimmed_fastqc.zip", sample=SAMPLES, read=["R1", "R2"])
    output:
        html = MULTIQC_DIR + "/multiqc_report.html"
    shell:
        """
        multiqc results/ -o {MULTIQC_DIR}
        """
