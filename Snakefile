# Import modules
import pandas as pd
import os

# CONFIG: Load config.yaml
configfile: "config/config.yaml"

# Load sample sheet
samples = pd.read_csv(config["samples"], sep="\t")
SAMPLES = samples["sample"].tolist()

# Directories from config
RAW_DIR = config["raw_dir"]
TRIMMED_DIR = config["trimmed_dir"]
FASTQC_DIR = config["fastqc_dir"]
TRIMMED_FASTQC_DIR = config["trimmed_fastqc_dir"]
MULTIQC_DIR = config["multiqc_dir"]
STAR_DIR = config["star_dir"]
COUNTS_DIR = config["counts_dir"]

# Define threads + parameters
STAR_INDEX = config["star_index"]
GTF = config["gtf"]
FASTA = config["fasta"]
STRAND = config["strand"]



# RULE: Define workflow
# define what output needs to be there to be done

rule all:
    input:
        # reference index
        os.path.join(STAR_INDEX, "SA"),
        # per-sample final BAM + index
        expand(os.path.join(STAR_DIR, "{sample}_Aligned.sortedByCoord.out.bam"), sample=SAMPLES),
        expand(os.path.join(STAR_DIR, "{sample}_Aligned.sortedByCoord.out.bam.bai"), sample=SAMPLES),
        # per-sample featureCounts
        expand(os.path.join(COUNTS_DIR, "{sample}_featureCounts.txt"), sample=SAMPLES),
        # combined count matrix
        os.path.join(COUNTS_DIR, "counts_combined.tsv"),
        # one MultiQC report for everything
        os.path.join(MULTIQC_DIR, "multiqc_report.html")

# QC
# RULE: FastQC on raw reads

rule fastqc_raw:
    input:
        fq = lambda wildcards: samples.loc[samples["sample"] == wildcards.sample, f"{wildcards.read}"].values[0]
    output:
        html = FASTQC_DIR + "/{sample}_{read}_fastqc.html",
        zip = temp(FASTQC_DIR + "/{sample}_{read}_fastqc.zip")
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
        R1_trimmed = temp(TRIMMED_DIR + "/{sample}_R1.trimmed.fastq.gz"),
        R2_trimmed = temp(TRIMMED_DIR + "/{sample}_R2.trimmed.fastq.gz"),
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
        zip = temp(TRIMMED_FASTQC_DIR + "/{sample}_{read}.trimmed_fastqc.zip")
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

# RULE: STAR indexing
# indexes the reference genome (once)

rule star_index:
  input:
    fasta = FASTA,
    gtf = GTF
  output:
    star_index = os.path.join(STAR_INDEX, "SA")
  threads: 8
  shell:
    """
    STAR --runThreadN {threads} \
         --runMode genomeGenerate \
         --genomeDir {STAR_INDEX} \
         --genomeFastaFiles {input.fasta} \
         --sjdbGTFfile {input.gtf} \
         --sjdbOverhang 150
    """


# RULE: STAR mapping
# align the filtered fastq files
# the fastq files from 1st and 2nd sequencing run are combined here to one BAM

rule star_align:
    input:
        R1 = os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fastq.gz"),
        R2 = os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fastq.gz"),
        star_index = os.path.join(STAR_INDEX, "SA")
    output:
        bam = os.path.join(STAR_DIR, "{sample}_Aligned.sortedByCoord.out.bam")
    threads: 8
    params:
        index = STAR_INDEX,
        gtf = GTF
    log:
        os.path.join(STAR_DIR, "{sample}_STAR.log")
    shell:
        """
        echo "PWD: $(pwd)"
        echo "Files: {input.R1} {input.R2}"
        ls -lh {input.R1}
        STAR \
          --runThreadN {threads} \
          --genomeDir {params.index} \
          --readFilesIn <(gunzip -c {input.R1}) <(gunzip -c {input.R2}) \
          --outFileNamePrefix {STAR_DIR}/{wildcards.sample}_ \
          --outSAMtype BAM SortedByCoordinate \
          --sjdbGTFfile {params.gtf} \
          --quantMode GeneCounts \
          > {log} 2>&1
        """

# RULE: index BAM file
# not necessarily needed for featureCounts 
# good for visualisation in genome browsers for example

rule index_bam:
    input:
        bam = os.path.join(STAR_DIR, "{sample}_Aligned.sortedByCoord.out.bam")
    output:
        bai = os.path.join(STAR_DIR, "{sample}_Aligned.sortedByCoord.out.bam.bai")
    shell:
        """
        samtools index {input.bam}
        """

# RULE: featureCounts
# quantify the mapped reads based on exons and sums up the reads for each gene
# featureCounts excludes reads that overlap multiple features 
# -O option would change this-one read can count for multiple genes

rule featureCounts:
    input:
        bam = os.path.join(STAR_DIR, "{sample}_Aligned.sortedByCoord.out.bam")
    output:
        counts = os.path.join(COUNTS_DIR, "{sample}_featureCounts.txt")
    threads: 4
    params:
        gtf = GTF,
        strand = STRAND
    log:
        os.path.join(COUNTS_DIR, "{sample}_featureCounts.log")
    shell:
        """
        featureCounts \
          -T {threads} \
          -a {params.gtf} \
          -o {output.counts} \
          -p -B -C -s {params.strand} \
          {input.bam} \
          > {log} 2>&1
        """

# RULE: create combined count matrix
# combine all reads in one count matrix with one column per sample

rule combine_counts:
    input:
        expand(os.path.join(COUNTS_DIR, "{sample}_featureCounts.txt"), sample=SAMPLES)
    output:
        combined = os.path.join(COUNTS_DIR, "counts_combined.tsv")
    run:
        import pandas as pd

        # Load all count files
        dfs = []
        for f in input:
            # Extract sample name from file name
            sample = os.path.basename(f).replace("_featureCounts.txt", "")
            # Read counts, skip the comment lines (start with '#')
            df = pd.read_csv(f, sep="\t", comment='#')
            # Keep only Geneid and counts column
            df = df[["Geneid", df.columns[-1]]]
            # Rename count column to sample name
            df = df.rename(columns={df.columns[-1]: sample})
            dfs.append(df)

        # Merge on Geneid
        combined = dfs[0]
        for df in dfs[1:]:
            combined = combined.merge(df, on="Geneid")

        # Save to output file
        combined.to_csv(output.combined, sep="\t", index=False)