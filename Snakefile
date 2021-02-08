__author__ = "Masood Zaka (https://github.com/masoodzaka/Snakemake_DNASeq.git)"
__licence__ = "MIT"

#shell.prefix("set -o pipefail; ")


configfile: "config.yaml"

SAMPLES, LANES = glob_wildcards("FASTQ/{sample}_{lane}_R1_001.fastq.gz")


rule all:
	input:
		"FastQC/multiqc_report.html",
		expand("SamToSortedBam/{sample}_{lane}.bam",zip, sample=SAMPLES, lane=LANES)

rule fastqc:
	input:
		R1=expand("FASTQ/{sample}_{lane}_R1_001.fastq.gz",zip, sample=SAMPLES, lane=LANES),
		R2=expand("FASTQ/{sample}_{lane}_R2_001.fastq.gz",zip, sample=SAMPLES, lane=LANES)

	output:
		expand("FastQC/{sample}_{lane}_R1_001_fastqc.html",zip, sample=SAMPLES, lane=LANES),
		expand("FastQC/{sample}_{lane}_R2_001_fastqc.html",zip, sample=SAMPLES, lane=LANES)

	params:
		threads = "5"

	shell:"""
		fastqc -t {params.threads} {input.R1} {input.R2} --noextract -q -o FastQC
	"""

rule multiqc:
	input:
		expand("FastQC/{sample}_{lane}_R1_001_fastqc.html",zip, sample=SAMPLES, lane=LANES),
		expand("FastQC/{sample}_{lane}_R2_001_fastqc.html",zip, sample=SAMPLES, lane=LANES)
	output:
		"FastQC/multiqc_report.html"
	shell:"""
			multiqc -f -q FastQC -o FastQC
		"""

rule bwa_mem:
	input:
		R1="FASTQ/{sample}_{lane}_R1_001.fastq.gz",
		R2="FASTQ/{sample}_{lane}_R2_001.fastq.gz"

	output:
		"BWA_MEM/{sample}_{lane}.sam"

	params:
		REF=config["REF"],
		RG=r"@RG\tID:{sample}\tLB:{sample}\tSM:{sample}\tPL:ILLUMINA"

	log:"LOGS/BWA_MEM/{sample}_{lane}.log"

	threads: 5

	message:"Running BWA_MEM for input file {input} using {threads} threads and saving as {output}"

	shell:"""
			bwa mem -t {threads} -M -v 1 -R '{params.RG}' {params.REF} {input.R1} {input.R2} > {output} 2> {log}
		"""

rule samtosortedbam:
	input:
		"BWA_MEM/{sample}_{lane}.sam"

	output:
		"SamToSortedBam/{sample}_{lane}.bam"

	params:
		REF=config["REF"],
		INDEX="true",
		VS="LENIENT",
		MRIR="3000000"

	log:"LOGS/SamToSortedBam/{sample}_{lane}.log"

	message:"Running SamSort for {input} saving sorted index BAM as {output}"

	shell:"""
			gatk --java-options "-Xmx14G -XX:ParallelGCThreads=2" SortSam \
			--INPUT {input} \
			--OUTPUT {output} \
			--CREATE_INDEX {params.INDEX} \
			-R {params.REF} \
			--VALIDATION_STRINGENCY {params.VS} \
			--MAX_RECORDS_IN_RAM {params.MRIR} \
			--SORT_ORDER coordinate 2> {log}
		"""

rule markduplicates:
	input:
		"SamToSortedBam/output.bam"

	output:
		"MarkDuplicates/output.dedup.bam"

	shell:"""
			gatk --java-options "-Xmx16G -XX:ParallelGCThreads=2"  MarkDuplicates \
			--INPUT {input} \
			--OUTPUT {output} \
			--METRICS_FILE MarkDuplicates/output.MD.metrics.txt \
			--CREATE_INDEX true \
			--VALIDATION_STRINGENCY STRICT \
			--MAX_RECORDS_IN_RAM 3000000
	"""

rule baserecal:
	input:
		"MarkDuplicates/output.dedup.bam"

	output:
		"BQSR_sample_lvl/output.Recal_data.grp"

	shell:"""
			gatk --java-options "-Xmx18G" BaseRecalibrator \
			-L intervals.interval_list \
			--interval-padding 100 \
			--input {input} \
			--known-sites /home/GATK_Bundle/b37/dbsnp_138.b37.vcf \
			--known-sites /home/GATK_Bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
			--known-sites /home/GATK_Bundle/b37/1000G_phase1.indels.b37.vcf \
			--reference /home/GATK_Bundle/b37/human_g1k_v37_decoy.fasta \
			--output {output}
	"""

rule applybqsr:
	input:
		"MarkDuplicates/output.dedup.bam"
	output:
		"BQSR_sample_lvl/output.recalibrated.bam"
	shell:"""
			gatk --java-options "-Xmx18G" ApplyBQSR \
			-L intervals.interval_list \
			--interval-padding 100 \
			-I MarkDuplicates/output.dedup.bam \
			-R /home/GATK_Bundle/b37/human_g1k_v37_decoy.fasta \
			--bqsr-recal-file BQSR_sample_lvl/output.Recal_data.grp \
			--output {output}
	"""