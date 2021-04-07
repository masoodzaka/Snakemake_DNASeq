__author__ = "Masood Zaka (https://github.com/masoodzaka/Snakemake_DNASeq.git)"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")

configfile: "config.yaml"


# include functions.py modules
include: "SCRIPTS/functions.py"

# Sample list
SAMPLES = MASTER_LIST["sample"].unique()

# Mutect list

NORMAL= list(MUTECT_LIST["normal"])
TUMOUR= list(MUTECT_LIST["tumour"])

# Target rules

ALL_RECALIBRATED_BAM = expand("BQSR_sample_lvl/{sample}.dedup.recalibrated.bam", sample=SAMPLES)
ALL_BAMFASTQC = expand("QC/BAMQC/{sample}.dedup.recalibrated_fastqc.html", sample=SAMPLES)
ALL_SAMTOOLSFLAGSTAT = expand("QC/SAMTOOLSFLAGSTAT/{sample}.dedup.recalibrated.flagstat", sample=SAMPLES)
ALL_HS_METRICS = expand("QC/HsMetrics/{sample}.dedup.recalibrated.hs_metrics.txt", sample=SAMPLES)
BAMQC = ["QC/BAMmultiqc_report.html"]
MUTECT = expand("MT2_Filt/{normal}.vs.{tumour}.somatic.vcf", zip,normal=NORMAL,tumour=TUMOUR)

ALL = []

ALL.extend(ALL_RECALIBRATED_BAM)
ALL.extend(ALL_BAMFASTQC)
ALL.extend(ALL_SAMTOOLSFLAGSTAT)
ALL.extend(ALL_HS_METRICS)
ALL.extend(BAMQC)
ALL.extend(MUTECT)

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
			bwa mem -t {threads} -M -v 1 {params.RG} {params.REF} {input.R1} {input.R2} > {output.SAM} 2> {log}
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
rule mutect_paired:
	input:
		unpack(mutect_inputs)

	output:
		BAM="MT2/{normal}.vs.{tumour}.bam",
		BAI="MT2/{normal}.vs.{tumour}.bai",
		F1R2="MT2/{normal}.vs.{tumour}.f1r2.tar.gz",
		VCF="MT2/{normal}.vs.{tumour}.vcf",
		IDX="MT2/{normal}.vs.{tumour}.vcf.idx",

	params:
		REF=config["REF"],
		INTERVALS=config["INTERVALS"],
		PADDING=config["PADDING"],
		AF_ONLY_GNOMAD=config["AF_ONLY_GNOMAD"],
		PON=config["PON"],
		B_NAME_N="{normal}"

	threads: 4

	log:"LOGS/MT2/{normal}.vs.{tumour}.log"

	benchmark:"LOGS/MT2/{normal}.vs.{tumour}.tsv"

	message: "Running GATK mutect2 for {input} using {threads} threads and saving as {output.BAM}"

	shell:"""
			gatk --java-options "-Xmx18G" Mutect2 \
			{params.INTERVALS} \
			--interval-padding {params.PADDING} \
			-I {input.TUMOUR} \
			-I {input.NORMAL} \
			-normal {params.B_NAME_N} \
			-R {params.REF} \
			--germline-resource {params.AF_ONLY_GNOMAD} \
			-pon {params.PON} \
			--f1r2-tar-gz {output.F1R2}\
			--output {output.VCF} \
			--bam-output {output.BAM} 2> {log}
	"""
rule readsOrientation:
	input:
		F1R2="MT2/{normal}.vs.{tumour}.f1r2.tar.gz",

	output:
		LROM="MT2_Filt/{normal}.vs.{tumour}.read-orientation-model.tar.gz",

	threads: 2

	log:"LOGS/MT2_Filt/{normal}.vs.{tumour}.readsOrientation.log"

	benchmark:"LOGS/MT2_Filt/{normal}.vs.{tumour}.readsOrientation.tsv"

	message: "Running GATK readsOrientation for {input} using {threads} threads and saving as {output}"

	shell:"""
			gatk --java-options "-Xmx2G" LearnReadOrientationModel \
			-I {input.F1R2} \
			-O {output.LROM} 2> {log}
	"""
rule getPileupsummaries:
	input:
		unpack(mutect_inputs),

	output:
		PS_N="MT2_Filt/{normal}.vs.{tumour}_N.getpileupsummaries.table",
		PS_T="MT2_Filt/{normal}.vs.{tumour}_T.getpileupsummaries.table",

	params:
		AF_ONLY_GNOMAD=config["AF_ONLY_GNOMAD"],

	threads: 2

	log:"LOGS/MT2_Filt/{normal}.vs.{tumour}.getPileupsummaries.log"

	benchmark:"LOGS/MT2_Filt/{normal}.vs.{tumour}.getPileupsummaries.tsv"

	message: "Running GATK getPileupsummaries for {input} using {threads} threads and saving as {output}"

	shell:"""
			gatk --java-options "-Xmx12G" GetPileupSummaries \
			-I {input.NORMAL} \
			-V {params.AF_ONLY_GNOMAD} \
			-L {params.AF_ONLY_GNOMAD} \
			-O {output.PS_N} 2> {log}

			gatk --java-options "-Xmx12G" GetPileupSummaries \
			-I {input.TUMOUR} \
			-V {params.AF_ONLY_GNOMAD} \
			-L {params.AF_ONLY_GNOMAD} \
			-O {output.PS_T} 2>> {log}

	"""
rule calculateContamination:
	input:
		PS_N="MT2_Filt/{normal}.vs.{tumour}_N.getpileupsummaries.table",
		PS_T="MT2_Filt/{normal}.vs.{tumour}_T.getpileupsummaries.table",

	output:
		CT="MT2_Filt/{normal}.vs.{tumour}.calculatecontamination.table",
		ST="MT2_Filt/{normal}.vs.{tumour}.segment.table",

	threads: 2

	log:"LOGS/MT2_Filt/{normal}.vs.{tumour}.calculateContamination.log"

	benchmark:"LOGS/MT2_Filt/{normal}.vs.{tumour}.calculateContamination.tsv"

	message: "Running GATK calculateContamination for {input} using {threads} threads and saving as {output}"

	shell:"""
			gatk --java-options "-Xmx12G" CalculateContamination \
			-I {input.PS_T}.getpileupsummaries.table \
			-matched {input.PS_N}.getpileupsummaries.table \
			-O {output.CT} \
			--tumor-segmentation {output.ST} 2> {log}
	"""
rule filter_mutectCalls_paired:
	input:
		F1R2="MT2/{normal}.vs.{tumour}.f1r2.tar.gz",
		VCF="MT2/{normal}.vs.{tumour}.vcf",
		IDX="MT2/{normal}.vs.{tumour}.vcf.idx",
		STATS="MT2/{normal}.vs.{tumour}.vcf.stats",
		CT="MT2_Filt/{normal}.vs.{tumour}.calculatecontamination.table",
		ST="MT2_Filt/{normal}.vs.{tumour}.segment.table",
		LROM="MT2_Filt/{normal}.vs.{tumour}.read-orientation-model.tar.gz",

	output:
		F_STATS="MT2_Filt/{normal}.vs.{tumour}.filtering.stats",
		Unfiltered_VCF="MT2_Filt/{normal}.vs.{tumour}.unfiltered.vcf",
		Unfiltered_IDX="MT2_Filt/{normal}.vs.{tumour}.unfiltered.vcf.idx",
		Somatic_VCF="MT2_Filt/{normal}.vs.{tumour}.somatic.vcf",
		Somatic_IDX="MT2_Filt/{normal}.vs.{tumour}.somatic.vcf.idx",

	params:
		REF=config["REF"],
		INTERVALS=config["INTERVALS"],
		PADDING=config["PADDING"],
		AF_ONLY_GNOMAD=config["AF_ONLY_GNOMAD"],
		B_NAME_N="{normal}",
		B_NAME_T="{tumour}"

	threads: 2

	log:"LOGS/MT2_Filt/{normal}.vs.{tumour}.log"

	benchmark:"LOGS/MT2_Filt/{normal}.vs.{tumour}.tsv"

	message: "Running GATK filtermutectCalls for {input.VCF} using {threads} threads and saving as {output.Somatic_VCF}"

	shell:"""
			gatk --java-options "-Xmx12G" FilterMutectCalls \
			-R {params.REF} \
			-V {input.VCF} \
			--tumor-segmentation {input.ST} \
			--contamination-table {input.CT} \
			--ob-priors {input.LROM} \
			-stats {input.STATS} \
			--filtering-stats {output.F_STATS} \
			-O {output.Unfiltered_VCF} 2> {log}

			gatk --java-options "-Xmx12G" SelectVariants \
			{params.INTERVALS} \
			--interval-padding {params.PADDING} \
			-V {output.Unfiltered_VCF} \
			-R {params.REF} \
			--exclude-filtered \
			--output {output.Somatic_VCF} 2>> {log}
	"""