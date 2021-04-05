__author__ = "Masood Zaka (https://github.com/masoodzaka/Snakemake_DNASeq.git)"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")

configfile: "config.yaml"

# include functions.py modules
include: "SCRIPTS/functions.py"

# Sample list
SAMPLES = MASTER_LIST["sample"].unique()

# Target rules

ALL_RECALIBRATED_BAM = expand("BQSR_sample_lvl/{sample}.dedup.recalibrated.bam", sample=SAMPLES)
ALL_BAMFASTQC = expand("QC/BAMQC/{sample}.dedup.recalibrated_fastqc.html", sample=SAMPLES)
ALL_SAMTOOLSFLAGSTAT = expand("QC/SAMTOOLSFLAGSTAT/{sample}.dedup.recalibrated.flagstat", sample=SAMPLES)
ALL_HS_METRICS = expand("QC/HsMetrics/{sample}.dedup.recalibrated.hs_metrics.txt", sample=SAMPLES)
BAMQC = ["QC/BAMmultiqc_report.html"]

ALL = []

ALL.extend(ALL_RECALIBRATED_BAM)
ALL.extend(ALL_BAMFASTQC)
ALL.extend(ALL_SAMTOOLSFLAGSTAT)
ALL.extend(ALL_HS_METRICS)
ALL.extend(BAMQC)

rule all:
	input:
		ALL

rule bwa_mem:
	input:
		unpack(bwa_input)
	output:
		SAM="BWA_MEM/{sample}.{runID}.sam"

	params:
		REF=config["REF"],
		RG=bwa_readgroup

	log:"LOGS/BWA_MEM/{sample}.{runID}.log"

	benchmark:"LOGS/BWA_MEM/{sample}.{runID}.tsv"

	threads: 5

	message:"Running BWA_MEM for input file {input} using {threads} threads and saving as {output}"

	shell:"""
			bwa mem -t {threads} -M -v 1 -R r'{params.RG}' {params.REF} {input.R1} {input.R2} > {output.SAM} 2> {log}
		"""

rule samtosortedbam:
	input:
		"BWA_MEM/{sample}.{runID}.sam"

	output:
		BAM="SamToSortedBam/{sample}.{runID}.bam"

	params:
		REF=config["REF"],
		INDEX="true",
		VS="LENIENT",
		MRIR="3000000",
		SO="coordinate"

	log:"LOGS/SamToSortedBam/{sample}.{runID}.log"

	benchmark:"LOGS/SamToSortedBam/{sample}.{runID}.tsv"

	message:"Running SamSort for {input} saving sorted index BAM as {output}"

	shell:"""
			gatk --java-options "-Xmx20G -XX:ParallelGCThreads=2" SortSam \
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
		mergebam_input,
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
			INPUT=""
			for i in {input}; do
				INPUT+=" --INPUT $i"
			done

			gatk --java-options "-Xmx20G -XX:ParallelGCThreads=2" MergeSamFiles \
			$INPUT \
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
			gatk --java-options "-Xmx20G -XX:ParallelGCThreads=2"  MarkDuplicates \
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

	message: "Running GATK BaseRecalibrator for {input} using {threads} threads and saving as {output}"

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
		BAM="BQSR_sample_lvl/{sample}.dedup.recalibrated.bam"

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

	message: "Running GATK ApplyBQSR for {input.BAM} using {threads} threads and saving as {output.BAM}"

	shell:"""
			gatk --java-options "-Xmx18G" ApplyBQSR \
			{params.INTERVALS} \
			--interval-padding {params.PADDING} \
			--input {input.BAM} \
			-R {params.REF} \
			--bqsr-recal-file {input.RECAL_DATA} \
			--output {output.BAM} 2> {log}
	"""
rule bamqc:
	input:
		"BQSR_sample_lvl/{sample}.dedup.recalibrated.bam"
	output:
		"QC/BAMQC/{sample}.dedup.recalibrated_fastqc.html",
		"QC/BAMQC/{sample}.dedup.recalibrated_fastqc.zip"

	params:
		BAM="bam"

	threads: 2

	log: "LOGS/QC/BAMQC/{sample}.log"

	benchmark: "LOGS/QC/BAMQC/{sample}.tsv"

	message: "Running BamQC for {input} and saving as {output}"

	shell:"""
		fastqc -f {params.BAM} -t {threads} --noextract {input} -o QC/BAMQC 2> {log}
	"""
rule samtoolsflagstat:
	input:
		"BQSR_sample_lvl/{sample}.dedup.recalibrated.bam"
	output:
		"QC/SAMTOOLSFLAGSTAT/{sample}.dedup.recalibrated.flagstat"

	threads: 1

	log: "LOGS/QC/SAMTOOLSFLAGSTAT/{sample}.log"

	benchmark: "LOGS/QC/SAMTOOLSFLAGSTAT/{sample}.tsv"

	message: "Running SamtoolsFlagstat for {input} and saving as {output}"

	shell:"""
		samtools flagstat {input} > {output} 2> {log}

	"""

rule HsMetrics:
	input:
		BAM="BQSR_sample_lvl/{sample}.dedup.recalibrated.bam"

	output:
		"QC/HsMetrics/{sample}.dedup.recalibrated.hs_metrics.txt"

	params:
		T_INTERVALS=config["T_INTERVALS"],
		B_INTERVALS=config["B_INTERVALS"]

	log: "LOGS/QC/HsMetrics/{sample}.HsMetrics.txt"

	benchmark: "LOGS/QC/HsMetrics/{sample}.HsMetrics.tsv"

	message: "Running GATK HsMetrics for {input} and saving as {output}"

	shell:"""
		gatk --java-options "-Xmx18G" CollectHsMetrics \
		--BAIT_INTERVALS {params.B_INTERVALS} \
		--INPUT {input.BAM} \
		--OUTPUT {output} \
		--TARGET_INTERVALS {params.T_INTERVALS} 2> {log}
	"""
rule multiqcBAM:
	input:
		multiqcbam_input

	output:
		"QC/BAMmultiqc_report.html"

	message: "Running MultiQCBAM for {input} and saving as {output}"

	shell:"""
		multiqc QC/SAMTOOLSFLAGSTAT QC/BAMQC QC/HsMetrics -n BAMmultiqc_report -d -f -q -o QC
	"""
