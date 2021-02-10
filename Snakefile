__author__ = "Masood Zaka (https://github.com/masoodzaka/Snakemake_DNASeq.git)"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")


configfile: "config.yaml"

SAMPLES, LANES = glob_wildcards("FASTQ/{sample}_{lane}_R1_001.fastq.gz")


rule all:
	input:
		"FastQC/multiqc_report.html",
		expand("SamToSortedBam/{sample}_{lane}.bam",zip, sample=SAMPLES, lane=LANES),
		expand("MergeBAMList/{sample}.bam", sample=SAMPLES),
		expand("MarkDuplicates/{sample}.dedup.bam", sample=SAMPLES),
		expand("BQSR_sample_lvl/{sample}.recalibrated.bam", sample=SAMPLES)

rule fastqc:
	input:
		R1=expand("FASTQ/{sample}_{lane}_R1_001.fastq.gz",zip, sample=SAMPLES, lane=LANES),
		R2=expand("FASTQ/{sample}_{lane}_R2_001.fastq.gz",zip, sample=SAMPLES, lane=LANES)

	output:
		expand("FastQC/{sample}_{lane}_R1_001_fastqc.html",zip, sample=SAMPLES, lane=LANES),
		expand("FastQC/{sample}_{lane}_R2_001_fastqc.html",zip, sample=SAMPLES, lane=LANES)

	threads: 5

	message:"Running FastQC for input file {input} using {threads} threads and saving as {output}"

	shell:"""
		fastqc -t {threads} {input.R1} {input.R2} --noextract -q -o FastQC
	"""

rule multiqc:
	input:
		expand("FastQC/{sample}_{lane}_R1_001_fastqc.html",zip, sample=SAMPLES, lane=LANES),
		expand("FastQC/{sample}_{lane}_R2_001_fastqc.html",zip, sample=SAMPLES, lane=LANES)

	output:
		"FastQC/multiqc_report.html"

	message:"Running MultiQC for input file {input} and saving as {output}"

	shell:"""
			multiqc -f -q FastQC -o FastQC
		"""

rule bwa_mem:
	input:
		R1="FASTQ/{sample}_{lane}_R1_001.fastq.gz",
		R2="FASTQ/{sample}_{lane}_R2_001.fastq.gz"

	output:
		SAM="BWA_MEM/{sample}_{lane}.sam"

	params:
		REF=config["REF"],
		RG=r"@RG\tID:{sample}\tLB:{sample}_{lane}\tSM:{sample}\tPL:ILLUMINA"

	log:"LOGS/BWA_MEM/{sample}_{lane}.log"

	benchmark:"LOGS/BWA_MEM/{sample}_{lane}.tsv"

	threads: 5

	message:"Running BWA_MEM for input file {input} using {threads} threads and saving as {output}"

	shell:"""
			bwa mem -t {threads} -M -v 1 -R '{params.RG}' {params.REF} {input.R1} {input.R2} > {output.SAM} 2> {log}
		"""

rule samtosortedbam:
	input:
		"BWA_MEM/{sample}_{lane}.sam"

	output:
		BAM="SamToSortedBam/{sample}_{lane}.bam"

	params:
		REF=config["REF"],
		INDEX="true",
		VS="LENIENT",
		MRIR="3000000",
		SO="coordinate"

	log:"LOGS/SamToSortedBam/{sample}_{lane}.log"

	benchmark:"LOGS/SamToSortedBam/{sample}_{lane}.tsv"

	message:"Running SamSort for {input} saving sorted index BAM as {output}"

	shell:"""
			gatk --java-options "-Xmx14G -XX:ParallelGCThreads=2" SortSam \
			--INPUT {input} \
			--OUTPUT {output.BAM} \
			--CREATE_INDEX {params.INDEX} \
			-R {params.REF} \
			--VALIDATION_STRINGENCY {params.VS} \
			--MAX_RECORDS_IN_RAM {params.MRIR} \
			--SORT_ORDER {params.SO} 2> {log}
		"""
rule mergebamlist:
	input:
		L1="SamToSortedBam/{sample}_L001.bam",
		L2="SamToSortedBam/{sample}_L002.bam"

	output:
		BAM="MergeBAMList/{sample}.bam"

	params:
		REF=config["REF"],
		INDEX="true",
		VS="STRICT",
		MRIR="8000000",
		THREADING="true",
		SO="coordinate"

	log:"LOGS/MergeBAMList/{sample}.log"

	benchmark:"LOGS/MergeBAMList/{sample}.tsv"

	message:"Running MergeSam for {input} saving as {output}"

	shell:"""
			gatk --java-options "-Xmx14G -XX:ParallelGCThreads=2" MergeSamFiles \
			--INPUT {input.L1} \
			--INPUT {input.L2} \
			--OUTPUT {output.BAM} \
			-R {params.REF} \
			--MAX_RECORDS_IN_RAM {params.MRIR} \
			--USE_THREADING {params.THREADING} \
			--SORT_ORDER {params.SO} \
			--CREATE_INDEX {params.INDEX} \
			--VALIDATION_STRINGENCY {params.VS} 2> {log}
	"""
rule markduplicates:
	input:
		"MergeBAMList/{sample}.bam"

	output:
		BAM="MarkDuplicates/{sample}.dedup.bam",
		METRICS="MarkDuplicates/{sample}.MD.metrics.txt"

	params:
		REF=config["REF"],
		INDEX="true",
		VS="STRICT",
		MRIR="3000000"

	log: "LOGS/MarkDuplicates/{sample}.log"

	benchmark:"LOGS/MarkDuplicates/{sample}.tsv"

	message:"Running MarkDuplicates for {input} saving as {output}"

	shell:"""
			gatk --java-options "-Xmx16G -XX:ParallelGCThreads=2"  MarkDuplicates \
			--INPUT {input} \
			--OUTPUT {output.BAM} \
			-R {params.REF} \
			--METRICS_FILE {output.METRICS} \
			--CREATE_INDEX {params.INDEX} \
			--VALIDATION_STRINGENCY {params.VS} \
			--MAX_RECORDS_IN_RAM {params.MRIR} 2> {log}
	"""

rule baserecal:
	input:
		"MarkDuplicates/{sample}.dedup.bam"

	output:
		RECAL_DATA="BQSR_sample_lvl/{sample}.Recal_data.grp"

	params:
		REF=config["REF"],
		INTERVALS=config["INTERVALS"],
		PADDING=config["PADDING"],
		DBSNP=config["DBSNP"],
		MILLS_1KG_GOLD=config["MILLS_1KG_GOLD"],
		PHASE1_INDELS=config["PHASE1_INDELS"],
		PHASE1_SNPS=config["PHASE1_SNPS"]

	threads: 4

	log:"LOGS/BQSR_sample_lvl/{sample}.recal.log"

	benchmark:"LOGS/BQSR_sample_lvl/{sample}.recal.tsv"

	message: "Running BaseRecalibrator for {input} using {threads} threads and saving as {output}"

	shell:"""
			gatk --java-options "-Xmx18G" BaseRecalibrator \
			-R {params.REF} \
			{params.INTERVALS} \
			--interval-padding {params.PADDING} \
			--input {input} \
			--known-sites {params.DBSNP} \
			--known-sites {params.MILLS_1KG_GOLD} \
			--known-sites {params.PHASE1_INDELS} \
			--known-sites {params.PHASE1_SNPS} \
			--output {output.RECAL_DATA} 2> {log}
	"""

rule applybqsr:
	input:
		BAM="MarkDuplicates/{sample}.dedup.bam",
		RECAL_DATA="BQSR_sample_lvl/{sample}.Recal_data.grp"

	output:
		BAM="BQSR_sample_lvl/{sample}.recalibrated.bam"

	params:
		REF=config["REF"],
		INTERVALS=config["INTERVALS"],
		PADDING=config["PADDING"],
		DBSNP=config["DBSNP"],
		MILLS_1KG_GOLD=config["MILLS_1KG_GOLD"],
		PHASE1_INDELS=config["PHASE1_INDELS"],
		PHASE1_SNPS=config["PHASE1_SNPS"]

	threads: 4

	log:"LOGS/BQSR_sample_lvl/{sample}.bqsr.log"

	benchmark:"LOGS/BQSR_sample_lvl/{sample}.bqsr.tsv"

	message: "Running ApplyBQSR for {input.BAM} using {threads} threads and saving as {output.BAM}"

	shell:"""
			gatk --java-options "-Xmx18G" ApplyBQSR \
			{params.INTERVALS} \
			--interval-padding {params.PADDING} \
			--input {input.BAM} \
			-R {params.REF} \
			--bqsr-recal-file {input.RECAL_DATA} \
			--output {output.BAM} 2> {log}
	"""