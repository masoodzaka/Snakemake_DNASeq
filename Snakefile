__author__ = "Masood Zaka (https://github.com/masoodzaka/Snakemake_DNASeq.git)"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")

configfile: "config.yaml"


# include functions.py modules
include: "SCRIPTS/functions.py"


# Target rules
#ALL_RECALIBRATED_BAM = expand("BQSR_sample_lvl/{sample}.dedup.recalibrated.bam", sample=SAMPLES)
#BAMQC = ["QC/BAMmultiqc_report.html"]
MUTECT_PAIRED = expand("MT2_Filt/{normal}.vs.{tumour}.somatic.vcf", zip, normal=MT2_Paired["normal"], tumour=MT2_Paired["tumour"])
VEP_PAIRED = expand("VEP/Paired/{normal}.vs.{tumour}.tsv", zip, normal=MT2_Paired["normal"], tumour=MT2_Paired["tumour"])
MUTECT_TUMOURONLY = expand("MT2_TumourOnly_Filt/{pon}.vs.{tumour}.somatic.vcf", pon=["PON"], tumour=MT2_TumourOnly)
VEP_TUMOURONLY = expand("VEP/TumourOnly/{pon}.vs.{tumour}.tsv", pon=["PON"], tumour=MT2_TumourOnly)


# extend the ALL rule using python extend list function
ALL = []

#ALL.extend(ALL_RECALIBRATED_BAM)
#ALL.extend(BAMQC)
if config["MUTECT2"]["Paired"]:
	ALL.extend(MUTECT_PAIRED)
if config["MUTECT2"]["TumourOnly"]:
	ALL.extend(MUTECT_TUMOURONLY)
ALL.extend(VEP_PAIRED)
ALL.extend(VEP_TUMOURONLY)

rule ALL:
	input:
		ALL,


rule BWA_Mem:
	input:
		unpack(bwa_input)
	output:
		SAM=temp("BWA_MEM/{sample}.{runID}.sam")

	params:
		REF=config["REF"],
		RG=bwa_readgroup

	log:"LOGS/BWA_MEM/{sample}.{runID}.log"

	benchmark:"LOGS/BWA_MEM/{sample}.{runID}.tsv"

	threads: 5

	conda : "ENVS/bwa.yaml"

	message:"Running BWA_MEM for input file {input} using {threads} threads and saving as {output}"

	shell:"""
			bwa mem -t {threads} -M -v 1 {params.RG} {params.REF} {input.R1} {input.R2} > {output.SAM} 2> {log}
		"""

rule SamToSortedBam:
	input:
		"BWA_MEM/{sample}.{runID}.sam"

	output:
		BAM=temp("SamToSortedBam/{sample}.{runID}.bam"),
		BAI=temp("SamToSortedBam/{sample}.{runID}.bai")

	params:
		REF=config["REF"],
		INDEX="true",
		VS="LENIENT",
		MRIR="3000000",
		SO="coordinate"

	log:"LOGS/SamToSortedBam/{sample}.{runID}.log"

	benchmark:"LOGS/SamToSortedBam/{sample}.{runID}.tsv"

	threads: 2

	conda: "ENVS/gatk4.yaml"

	message:"Running SamSort for {input} using {threads} threads and saving sorted index BAM as {output}"

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

rule MergeBamList:
	input:
		mergebam_input,
	output:
		BAM=temp("MergeBAMList/{sample}.bam"),
		BAI=temp("MergeBAMList/{sample}.bai")

	params:
		REF=config["REF"],
		INDEX="true",
		VS="STRICT",
		MRIR="8000000",
		THREADING="true",
		SO="coordinate"

	log:"LOGS/MergeBAMList/{sample}.log"

	benchmark:"LOGS/MergeBAMList/{sample}.tsv"

	threads: 2

	conda: "ENVS/gatk4.yaml",

	message:"Running MergeSam for {input} using {threads} threads and using {threads} saving as {output}"

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
rule MarkDuplicates:
	input:
		BAM="MergeBAMList/{sample}.bam",
		BAI="MergeBAMList/{sample}.bai"

	output:
		BAM=temp("MarkDuplicates/{sample}.dedup.bam"),
		BAI=temp("MarkDuplicates/{sample}.dedup.bai"),
		METRICS="MarkDuplicates/{sample}.MD.metrics.txt"

	params:
		REF=config["REF"],
		INDEX="true",
		VS="STRICT",
		MRIR="3000000"

	log: "LOGS/MarkDuplicates/{sample}.log"

	benchmark:"LOGS/MarkDuplicates/{sample}.tsv"

	threads: 2

	conda: "ENVS/gatk4.yaml",

	message:"Running MarkDuplicates for {input} using {threads} threads and saving as {output}"

	shell:"""
			gatk --java-options "-Xmx20G -XX:ParallelGCThreads=2"  MarkDuplicates \
			--INPUT {input.BAM} \
			--OUTPUT {output.BAM} \
			-R {params.REF} \
			--METRICS_FILE {output.METRICS} \
			--CREATE_INDEX {params.INDEX} \
			--VALIDATION_STRINGENCY {params.VS} \
			--MAX_RECORDS_IN_RAM {params.MRIR} 2> {log}
	"""
rule BaseRecal:
	input:
		BAM="MarkDuplicates/{sample}.dedup.bam",
		BAI="MarkDuplicates/{sample}.dedup.bai",


	output:
		RECAL_DATA="BQSR_sample_lvl/{sample}.Recal_data.grp"

	params:
		REF=config["REF"],
		INTERVALS=config["INTERVALS"] if config["SEQUENCING"]["WES"] else " ",
		PADDING=config["PADDING"] if config["SEQUENCING"]["WES"] else " ",
		DBSNP=config["DBSNP"],
		MILLS_1KG_GOLD=config["MILLS_1KG_GOLD"],
		PHASE1_INDELS=config["PHASE1_INDELS"],
		PHASE1_SNPS=config["PHASE1_SNPS"]

	threads: 4

	conda: "ENVS/gatk4.yaml",

	log:"LOGS/BQSR_sample_lvl/{sample}.recal.log"

	benchmark:"LOGS/BQSR_sample_lvl/{sample}.recal.tsv"

	message: "Running GATK BaseRecalibrator for {input} using {threads} threads and saving as {output}"

	shell:"""
			gatk --java-options "-Xmx18G" BaseRecalibrator \
			-R {params.REF} \
			{params.INTERVALS} \
			--interval-padding {params.PADDING} \
			--input {input.BAM} \
			--known-sites {params.DBSNP} \
			--known-sites {params.MILLS_1KG_GOLD} \
			--known-sites {params.PHASE1_INDELS} \
			--known-sites {params.PHASE1_SNPS} \
			--output {output.RECAL_DATA} 2> {log}
	"""

rule ApplyBqsr:
	input:
		BAM="MarkDuplicates/{sample}.dedup.bam",
		BAI="MarkDuplicates/{sample}.dedup.bai",
		RECAL_DATA="BQSR_sample_lvl/{sample}.Recal_data.grp"

	output:
		BAM=protected("BQSR_sample_lvl/{sample}.dedup.recalibrated.bam"),
		BAI=protected("BQSR_sample_lvl/{sample}.dedup.recalibrated.bai")

	params:
		REF=config["REF"],
		INTERVALS= config["INTERVALS"] if config["SEQUENCING"]["WES"] else " ",
		PADDING=config["PADDING"] if config["SEQUENCING"]["WES"] else 0,
		DBSNP=config["DBSNP"],
		MILLS_1KG_GOLD=config["MILLS_1KG_GOLD"],
		PHASE1_INDELS=config["PHASE1_INDELS"],
		PHASE1_SNPS=config["PHASE1_SNPS"]

	threads: 4

	conda: "ENVS/gatk4.yaml",

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
rule BamQC:
	input:
		"BQSR_sample_lvl/{sample}.dedup.recalibrated.bam"
	output:
		"QC/BAMQC/{sample}.dedup.recalibrated_fastqc.html",
		"QC/BAMQC/{sample}.dedup.recalibrated_fastqc.zip"

	params:
		BAM="bam"

	threads: 2

	conda: "ENVS/qc.yaml",

	log: "LOGS/QC/BAMQC/{sample}.log"

	benchmark: "LOGS/QC/BAMQC/{sample}.tsv"

	message: "Running BamQC for {input} using {threads} threads and saving as {output}"

	shell:"""
		fastqc -f {params.BAM} -t {threads} --noextract {input} -o QC/BAMQC 2> {log}
	"""
rule SamtoolsFlagStat:
	input:
		"BQSR_sample_lvl/{sample}.dedup.recalibrated.bam"
	output:
		"QC/SAMTOOLSFLAGSTAT/{sample}.dedup.recalibrated.flagstat"

	threads: 1

	conda: "ENVS/samtools.yaml",

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
		T_INTERVALS=config["T_INTERVALS"] if config["SEQUENCING"]["WES"] else " ",
		B_INTERVALS=config["B_INTERVALS"] if config["SEQUENCING"]["WES"] else 0,

	log: "LOGS/QC/HsMetrics/{sample}.HsMetrics.txt"

	threads: 4,

	conda: "ENVS/gatk4.yaml",

	benchmark: "LOGS/QC/HsMetrics/{sample}.HsMetrics.tsv"

	message: "Running GATK HsMetrics for {input} using {threads} threads and saving as {output}"

	shell:"""
		gatk --java-options "-Xmx18G" CollectHsMetrics \
		--BAIT_INTERVALS {params.B_INTERVALS} \
		--INPUT {input.BAM} \
		--OUTPUT {output} \
		--TARGET_INTERVALS {params.T_INTERVALS} 2> {log}
	"""
rule MultiqcBAM:
	input:
		multiqcbam_input

	output:
		"QC/BAMmultiqc_report.html"

	message: "Running MultiQCBAM for {input} using {threads} threads and saving as {output}"

	threads: 2,

	conda: "ENVS/qc.yaml"

	shell:"""
		multiqc QC/SAMTOOLSFLAGSTAT QC/BAMQC QC/HsMetrics -n BAMmultiqc_report -d -f -q -o QC
	"""
rule Mutect_Paired:
	input:
		unpack(mutect_paired_inputs)

	output:
		BAM=protected("MT2/{normal}.vs.{tumour}.bam"),
		BAI=protected("MT2/{normal}.vs.{tumour}.bai"),
		F1R2=temp("MT2/{normal}.vs.{tumour}.f1r2.tar.gz"),
		VCF=temp("MT2/{normal}.vs.{tumour}.vcf"),
		IDX=temp("MT2/{normal}.vs.{tumour}.vcf.idx"),
		STATS=temp("MT2/{normal}.vs.{tumour}.vcf.stats"),

	params:
		REF=config["REF"],
		INTERVALS=config["INTERVALS"] if config["SEQUENCING"]["WES"] else " ",
		PADDING=config["PADDING"] if config["SEQUENCING"]["WES"] else 0,
		AF_ONLY_GNOMAD=config["AF_ONLY_GNOMAD"],
		PON=config["PON"],
		B_NAME_N="{normal}"

	threads: 4

	conda: "ENVS/gatk4.yaml",

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
rule ReadsOrientation_Paired:
	input:
		F1R2="MT2/{normal}.vs.{tumour}.f1r2.tar.gz",

	output:
		LROM=temp("MT2_Filt/{normal}.vs.{tumour}.read-orientation-model.tar.gz"),

	threads: 2

	conda: "ENVS/gatk4.yaml",

	log:"LOGS/MT2_Filt/{normal}.vs.{tumour}.readsOrientation.log"

	benchmark:"LOGS/MT2_Filt/{normal}.vs.{tumour}.readsOrientation.tsv"

	message: "Running GATK readsOrientation for {input} using {threads} threads and saving as {output}"

	shell:"""
			gatk --java-options "-Xmx2G" LearnReadOrientationModel \
			-I {input.F1R2} \
			-O {output.LROM} 2> {log}
	"""
rule GetPileupsummaries_Paired:
	input:
		unpack(mutect_paired_inputs),

	output:
		PS_N=temp("MT2_Filt/{normal}.vs.{tumour}_N.getpileupsummaries.table"),
		PS_T=temp("MT2_Filt/{normal}.vs.{tumour}_T.getpileupsummaries.table"),

	params:
		AF_ONLY_GNOMAD=config["AF_ONLY_GNOMAD"],

	threads: 2

	conda: "ENVS/gatk4.yaml",

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
rule CalculateContamination_Paired:
	input:
		PS_N="MT2_Filt/{normal}.vs.{tumour}_N.getpileupsummaries.table",
		PS_T="MT2_Filt/{normal}.vs.{tumour}_T.getpileupsummaries.table",

	output:
		CT=temp("MT2_Filt/{normal}.vs.{tumour}.calculatecontamination.table"),
		ST=temp("MT2_Filt/{normal}.vs.{tumour}.segment.table"),

	threads: 2

	conda: "ENVS/gatk4.yaml",

	log:"LOGS/MT2_Filt/{normal}.vs.{tumour}.calculateContamination.log"

	benchmark:"LOGS/MT2_Filt/{normal}.vs.{tumour}.calculateContamination.tsv"

	message: "Running GATK calculateContamination for {input} using {threads} threads and saving as {output}"

	shell:"""
			gatk --java-options "-Xmx12G" CalculateContamination \
			-I {input.PS_T} \
			-matched {input.PS_N} \
			-O {output.CT} \
			--tumor-segmentation {output.ST} 2> {log}
	"""
rule Filter_MutectCalls_Paired:
	input:
		F1R2="MT2/{normal}.vs.{tumour}.f1r2.tar.gz",
		VCF="MT2/{normal}.vs.{tumour}.vcf",
		IDX="MT2/{normal}.vs.{tumour}.vcf.idx",
		STATS="MT2/{normal}.vs.{tumour}.vcf.stats",
		CT="MT2_Filt/{normal}.vs.{tumour}.calculatecontamination.table",
		ST="MT2_Filt/{normal}.vs.{tumour}.segment.table",
		LROM="MT2_Filt/{normal}.vs.{tumour}.read-orientation-model.tar.gz",

	output:
		F_STATS=protected("MT2_Filt/{normal}.vs.{tumour}.filtering.stats"),
		Unfiltered_VCF=protected("MT2_Filt/{normal}.vs.{tumour}.unfiltered.vcf"),
		Unfiltered_IDX=protected("MT2_Filt/{normal}.vs.{tumour}.unfiltered.vcf.idx"),
		Somatic_VCF=protected("MT2_Filt/{normal}.vs.{tumour}.somatic.vcf"),
		Somatic_IDX=protected("MT2_Filt/{normal}.vs.{tumour}.somatic.vcf.idx"),

	params:
		REF=config["REF"],
		INTERVALS=config["INTERVALS"] if config["SEQUENCING"]["WES"] else " ",
		PADDING=config["PADDING"] if config["SEQUENCING"]["WES"] else 0,
		AF_ONLY_GNOMAD=config["AF_ONLY_GNOMAD"],
		B_NAME_N="{normal}",
		B_NAME_T="{tumour}"

	threads: 2

	conda: "ENVS/gatk4.yaml",

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

rule VEP_Paired:
	input:
		VCF="MT2_Filt/{normal}.vs.{tumour}.somatic.vcf",
		IDX="MT2_Filt/{normal}.vs.{tumour}.somatic.vcf.idx",

	output:
		VEP=protected("VEP/Paired/{normal}.vs.{tumour}.tsv"),

	params:
		REF=config["REF"],
		VEP_DATA=config["VEP_DATA"],

	threads: 2

	conda: "ENVS/vep.yaml"

	log:"LOGS/VEP/{normal}.vs.{tumour}.log"

	benchmark:"LOGS/VEP/{normal}.vs.{tumour}.tsv"

	message: "Running Ensembl VEP for {input.VCF} using {threads} threads and saving as {output.VEP}"

	shell:"""
			vep \
			-i {input.VCF} \
			--offline \
			--dir {params.VEP_DATA} \
			--dir_cache {params.VEP_DATA} \
			--everything \
			--nearest symbol \
			--total_length \
			--force_overwrite \
			--tab \
			--plugin ExACpLI,{params.VEP_DATA}/Plugins/ExACpLI_values.txt \
			--plugin LoFtool,{params.VEP_DATA}/Plugins/LoFtool_scores.txt \
			--plugin Carol \
			--plugin Blosum62 \
			-o {output.VEP} \
			--buffer_size 5000 \
			--fork 10 \
			--pick_allele \
			--show_ref_allele 2>{log}
	"""

rule Mutect_TumorOnly:
	input:
		unpack(mutect_tumourOnly_inputs)

	output:
		BAM=protected("MT2_TumourOnly/{pon}.vs.{tumour}.bam"),
		BAI=protected("MT2_TumourOnly/{pon}.vs.{tumour}.bai"),
		F1R2=temp("MT2_TumourOnly/{pon}.vs.{tumour}.f1r2.tar.gz"),
		VCF=temp("MT2_TumourOnly/{pon}.vs.{tumour}.vcf"),
		IDX=temp("MT2_TumourOnly/{pon}.vs.{tumour}.vcf.idx"),
		STATS=temp("MT2_TumourOnly/{pon}.vs.{tumour}.vcf.stats"),

	params:
		REF=config["REF"],
		INTERVALS=config["INTERVALS"] if config["SEQUENCING"]["WES"] else " ",
		PADDING=config["PADDING"] if config["SEQUENCING"]["WES"] else 0,
		AF_ONLY_GNOMAD=config["AF_ONLY_GNOMAD"],
		PON=config["PON"],

	threads: 4

	conda: "ENVS/gatk4.yaml",

	log:"LOGS/MT2_TumourOnly/{pon}.vs.{tumour}.log"

	benchmark:"LOGS/MT2_TumourOnly/{pon}.vs.{tumour}.tsv"

	message: "Running GATK mutect2 for {input} using {threads} threads and saving as {output.BAM}"

	shell:"""
			gatk --java-options "-Xmx18G" Mutect2 \
			{params.INTERVALS} \
			--interval-padding {params.PADDING} \
			-I {input.TUMOUR} \
			-R {params.REF} \
			--germline-resource {params.AF_ONLY_GNOMAD} \
			-pon {params.PON} \
			--f1r2-tar-gz {output.F1R2}\
			--output {output.VCF} \
			--bam-output {output.BAM} 2> {log}
	"""
rule ReadsOrientation_TumourOnly:
	input:
		F1R2="MT2_TumourOnly/{pon}.vs.{tumour}.f1r2.tar.gz",

	output:
		LROM=temp("MT2_TumourOnly_Filt/{pon}.vs.{tumour}.read-orientation-model.tar.gz"),

	threads: 2

	conda: "ENVS/gatk4.yaml",

	log:"LOGS/MT2_TumourOnly_Filt/{pon}.vs.{tumour}.readsOrientation_TumourOnly.log"

	benchmark:"LOGS/MT2_TumourOnly_Filt/{pon}.vs.{tumour}.readsOrientation_TumourOnly.tsv"

	message: "Running GATK readsOrientation_TumourOnly for {input} using {threads} threads and saving as {output}"

	shell:"""
			gatk --java-options "-Xmx2G" LearnReadOrientationModel \
			-I {input.F1R2} \
			-O {output.LROM} 2> {log}
	"""
rule GetPileupsummaries_TumourOnly:
	input:
		unpack(mutect_tumourOnly_inputs)

	output:
		PS_T=temp("MT2_TumourOnly_Filt/{pon}.vs.{tumour}_T.getpileupsummaries.table"),

	params:
		AF_ONLY_GNOMAD=config["AF_ONLY_GNOMAD"],

	threads: 2

	conda: "ENVS/gatk4.yaml",

	log:"LOGS/MT2_TumourOnly_Filt/{pon}.vs.{tumour}.getPileupsummaries_TumourOnly.log"

	benchmark:"LOGS/MT2_TumourOnly_Filt/{pon}.vs.{tumour}.getPileupsummaries_TumourOnly.tsv"

	message: "Running GATK getPileupsummaries_TumourOnly for {input} using {threads} threads and saving as {output}"

	shell:"""

			gatk --java-options "-Xmx12G" GetPileupSummaries \
			-I {input.TUMOUR} \
			-V {params.AF_ONLY_GNOMAD} \
			-L {params.AF_ONLY_GNOMAD} \
			-O {output.PS_T} 2>> {log}

	"""
rule CalculateContamination_TumourOnly:
	input:
		PS_T="MT2_TumourOnly_Filt/{pon}.vs.{tumour}_T.getpileupsummaries.table",

	output:
		CT=temp("MT2_TumourOnly_Filt/{pon}.vs.{tumour}.calculateContamination_TumourOnly.table"),
		ST=temp("MT2_TumourOnly_Filt/{pon}.vs.{tumour}.segment_TumourOnly.table"),

	threads: 2

	conda: "ENVS/gatk4.yaml",

	log:"LOGS/MT2_TumourOnly_Filt/{pon}.vs.{tumour}.calculateContamination_TumourOnly.log"

	benchmark:"LOGS/MT2_TumourOnly_Filt/{pon}.vs.{tumour}.calculateContamination_TumourOnly.tsv"

	message: "Running GATK calculateContamination_TumourOnly for {input} using {threads} threads and saving as {output}"

	shell:"""
			gatk --java-options "-Xmx12G" CalculateContamination \
			-I {input.PS_T} \
			-O {output.CT} \
			--tumor-segmentation {output.ST} 2> {log}
	"""
rule Filter_MutectCalls_TumourOly:
	input:
		F1R2="MT2_TumourOnly/{pon}.vs.{tumour}.f1r2.tar.gz",
		VCF="MT2_TumourOnly/{pon}.vs.{tumour}.vcf",
		IDX="MT2_TumourOnly/{pon}.vs.{tumour}.vcf.idx",
		STATS="MT2_TumourOnly/{pon}.vs.{tumour}.vcf.stats",
		CT="MT2_TumourOnly_Filt/{pon}.vs.{tumour}.calculateContamination_TumourOnly.table",
		ST="MT2_TumourOnly_Filt/{pon}.vs.{tumour}.segment_TumourOnly.table",
		LROM="MT2_TumourOnly_Filt/{pon}.vs.{tumour}.read-orientation-model.tar.gz",

	output:
		F_STATS=protected("MT2_TumourOnly_Filt/{pon}.vs.{tumour}.filtering.stats"),
		Unfiltered_VCF=protected("MT2_TumourOnly_Filt/{pon}.vs.{tumour}.unfiltered.vcf"),
		Unfiltered_IDX=protected("MT2_TumourOnly_Filt/{pon}.vs.{tumour}.unfiltered.vcf.idx"),
		Somatic_VCF=protected("MT2_TumourOnly_Filt/{pon}.vs.{tumour}.somatic.vcf"),
		Somatic_IDX=protected("MT2_TumourOnly_Filt/{pon}.vs.{tumour}.somatic.vcf.idx"),

	params:
		REF=config["REF"],
		INTERVALS=config["INTERVALS"] if config["SEQUENCING"]["WES"] else " ",
		PADDING=config["PADDING"] if config["SEQUENCING"]["WES"] else 0,
		AF_ONLY_GNOMAD=config["AF_ONLY_GNOMAD"],

	threads: 2

	conda: "ENVS/gatk4.yaml",

	log:"LOGS/MT2_TumourOnly_Filt/{pon}.vs.{tumour}.log"

	benchmark:"LOGS/MT2_TumourOnly_Filt/{pon}.vs.{tumour}.tsv"

	message: "Running GATK Filter_MutectCalls_TumourOly for {input.VCF} using {threads} threads and saving as {output.Somatic_VCF}"

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

rule VEP_TumourOnly:
	input:
		VCF="MT2_TumourOnly_Filt/{pon}.vs.{tumour}.somatic.vcf",
		IDX="MT2_TumourOnly_Filt/{pon}.vs.{tumour}.somatic.vcf.idx",

	output:
		VEP=protected("VEP/TumourOnly/{pon}.vs.{tumour}.tsv"),

	params:
		REF=config["REF"],
		VEP_DATA=config["VEP_DATA"],

	threads: 2

	conda: "ENVS/vep.yaml"

	log:"LOGS/VEP/{pon}.vs.{tumour}.log"

	benchmark:"LOGS/VEP/{pon}.vs.{tumour}.tsv"

	message: "Running Ensembl VEP for {input.VCF} using {threads} threads and saving as {output.VEP}"

	shell:"""
			vep \
			-i {input.VCF} \
			--offline \
			--dir {params.VEP_DATA} \
			--dir_cache {params.VEP_DATA} \
			--everything \
			--nearest symbol \
			--total_length \
			--force_overwrite \
			--tab \
			--plugin ExACpLI,{params.VEP_DATA}/Plugins/ExACpLI_values.txt \
			--plugin LoFtool,{params.VEP_DATA}/Plugins/LoFtool_scores.txt \
			--plugin Carol \
			--plugin Blosum62 \
			-o {output.VEP} \
			--buffer_size 5000 \
			--fork 10 \
			--pick_allele \
			--show_ref_allele 2>{log}
	"""
