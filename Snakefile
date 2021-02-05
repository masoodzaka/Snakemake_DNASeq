rule fastqc:
	input:
		"data/samples/S70_L001_R1_001.fastq.gz",
		"data/samples/S70_L001_R2_001.fastq.gz"
	output:
		"FastQC/S70_L001_R1_001_fastqc.html",
		"FastQC/S70_L001_R1_001_fastqc.zip",
		"FastQC/S70_L001_R2_001_fastqc.html",
		"FastQC/S70_L001_R2_001_fastqc.zip"
	shell:"""
			fastqc -t 4 --noextract -q data/samples/S70_L001_R1_001.fastq.gz data/samples/S70_L001_R2_001.fastq.gz -o FastQC
		"""
rule multiqc:
	input:
		"FastQC/S70_L001_R1_001_fastqc.html",
		"FastQC/S70_L001_R2_001_fastqc.html"
	output:
		"multiqc_report.html",
		directory("multiqc_data")
	shell:"""
			multiqc -f -q FastQC
		"""
rule bwa_mem:
	input:
		"/home/GATK_Bundle/b37/human_g1k_v37_decoy.fasta",
		"data/samples/S70_L001_R1_001.fastq.gz",
		"data/samples/S70_L001_R2_001.fastq.gz"
	output:
		"BWA_MEM/output.sam"
	params:
		rg=r"@RG\tID:RG_01\tLB:Lib_01\tSM:Sample_01\tPL:ILLUMINA"
	shell:"""
			bwa mem -t 5 -M -v 1 -R '{params.rg}' {input} > {output}
		"""
rule samtosortedbam:
	shell:"""
			gatk --java-options "-Xmx14G -XX:ParallelGCThreads=2" SortSam \
			--INPUT BWA_MEM/output.sam \
			--OUTPUT SamToSortedBam/output.bam \
			--CREATE_INDEX true \
			--VALIDATION_STRINGENCY LENIENT \
			--MAX_RECORDS_IN_RAM 3000000 \
			--SORT_ORDER coordinate
		"""
rule markduplicates:
	shell:"""
			gatk --java-options "-Xmx16G -XX:ParallelGCThreads=2"  MarkDuplicates \
			--INPUT SamToSortedBam/output.bam \
			--OUTPUT MarkDuplicates/output.dedup.bam \
			--METRICS_FILE MarkDuplicates/output.MD.metrics.txt \
			--CREATE_INDEX true \
			--VALIDATION_STRINGENCY STRICT \
			--MAX_RECORDS_IN_RAM 3000000
	"""
rule baserecal:
	shell:"""
			gatk --java-options "-Xmx18G" BaseRecalibrator \
			-L intervals.interval_list \
			--interval-padding 100 \
			--input MarkDuplicates/output.dedup.bam \
			--known-sites /home/GATK_Bundle/b37/dbsnp_138.b37.vcf \
			--known-sites /home/GATK_Bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
			--known-sites /home/GATK_Bundle/b37/1000G_phase1.indels.b37.vcf \
			--reference /home/GATK_Bundle/b37/human_g1k_v37_decoy.fasta \
			--output BQSR_sample_lvl/output.Recal_data.grp
	"""
rule applybqsr:
	output:
			"output.recalibrated.bam"
	shell:"""
			gatk --java-options "-Xmx18G" ApplyBQSR \
			-L intervals.interval_list \
			--interval-padding 100 \
			-I MarkDuplicates/output.dedup.bam \
			-R /home/GATK_Bundle/b37/human_g1k_v37_decoy.fasta \
			--bqsr-recal-file BQSR_sample_level/output.Recal_data.grp \
			--output BQSR_sample_level/output.recalibrated.bam
	"""