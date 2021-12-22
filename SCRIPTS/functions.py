__author__ = "Masood Zaka (https://github.com/masoodzaka/Snakemake_DNASeq.git)"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

configfile: "config.yaml"

MASTER_LIST = pd.read_table(config["MASTER_LIST"], dtype=str).set_index(["sample","runID"], drop=False)
MASTER_LIST.index = MASTER_LIST.index.set_levels([i.astype(str) for i in MASTER_LIST.index.levels])

MUTECT_LIST = pd.read_table(config["MUTECT_LIST"]).set_index(["normal","tumour"], drop=False)
MUTECT_LIST.index = MUTECT_LIST.index.set_levels([i.astype(str) for i in MUTECT_LIST.index.levels])

# samples list from master_list

SAMPLES = [ sample for sample in MASTER_LIST["sample"]]

# drop nan from from mutect list to get the paired samples only
MT2_Paired = (MUTECT_LIST).dropna()

# aggregate tumour samples using tumour column of mutect list
MT2_TumourOnly = [tumour for tumour in MUTECT_LIST["tumour"]]

# functions for rules inputs

def bwa_input(wildcards):
	""" function provides the rowwise fastqs from the master list in the form a python dictionary"""
	fastqs = MASTER_LIST.loc[(wildcards.sample, wildcards.runID), ["fastq1", "fastq2"]].dropna()
	if len(fastqs) == 2:
		return {"R1": fastqs.fastq1, "R2": fastqs.fastq2}
	else:
		return {"R1": fastqs.fastq1}

def bwa_readgroup(wildcards):
	"""this functions return the read group column for each sample and lane run"""
	return r"-R '{rg.rg}'".format(rg=MASTER_LIST.loc[(wildcards.sample, wildcards.runID), ["rg"]].dropna())


def mergebam_input(wildcards):
	return expand("SamToSortedBam/{sample}.{runID}.bam",
		sample=wildcards.sample,
		runID=MASTER_LIST.loc[wildcards.sample].runID)

def multiqcbam_input(wildcards):
	return expand(["QC/BAMQC/{sample}.dedup.recalibrated_fastqc.zip",
		"QC/SAMTOOLSFLAGSTAT/{sample}.dedup.recalibrated.flagstat",
		"QC/HsMetrics/{sample}.dedup.recalibrated.hs_metrics.txt"],
		sample=SAMPLES)

def mutect_paired_inputs(wildcards):
	inputs = MT2_Paired.loc[(wildcards.normal, wildcards.tumour), ["normal", "tumour"]].dropna()
	return {"NORMAL": "BQSR_sample_lvl/{}.dedup.recalibrated.bam".format(inputs.normal),"TUMOUR": "BQSR_sample_lvl/{}.dedup.recalibrated.bam".format(inputs.tumour)}

def mutect_tumourOnly_inputs(wildcards):
	return {"TUMOUR": "BQSR_sample_lvl/{}.dedup.recalibrated.bam".format(tumour) for tumour in MT2_TumourOnly}