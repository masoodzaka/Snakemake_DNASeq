__author__ = "Masood Zaka (https://github.com/masoodzaka/Snakemake_DNASeq.git)"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

configfile: "config.yaml"

MASTER_LIST = pd.read_table(config["MASTER_LIST"], dtype=str).set_index(["sample","runID"], drop=False)
MASTER_LIST.index = MASTER_LIST.index.set_levels([i.astype(str) for i in MASTER_LIST.index.levels])

##### Wildcard constraints #####
wildcard_constraints:
	vartype="snvs|indels",
	sample="|".join(MASTER_LIST["sample"].unique()),
	runID="|".join(MASTER_LIST["runID"])

def bwa_input(wildcards):
	fastqs = MASTER_LIST.loc[(wildcards.sample, wildcards.runID), ["fastq1", "fastq2"]].dropna()
	if len(fastqs) == 2:
		return {"R1": fastqs.fastq1, "R2": fastqs.fastq2}
	return {"R1": fastqs.fastq1}

def bwa_readgroup(wildcards):
	rg = MASTER_LIST.loc[(wildcards.sample, wildcards.runID), ["rg"]].dropna()
	return rg.rg

def mergebam_input(wildcards):
	return expand("SamToSortedBam/{sample}.{runID}.bam",
		sample=wildcards.sample,
		runID=MASTER_LIST.loc[wildcards.sample].runID)

def multiqcbam_input(wildcards):
	return expand(["QC/SAMTOOLSFLAGSTAT/{sample}.dedup.recalibrated.flagstat",
		"QC/BAMQC/{sample}.dedup.recalibrated_fastqc.zip",
		"QC/HsMetrics/{sample}.dedup.recalibrated.hs_metrics.txt"],
		sample=wildcards.sample)