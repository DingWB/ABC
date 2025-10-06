import os
import pandas as pd
import hashlib
import numpy as np
from snakemake.utils import min_version
from functools import partial
# min_version("7.0")
import time
import random
from pandas.errors import EmptyDataError

configfile: "config/config.yaml"
SNAKEFILE_DIR = workflow.basedir

# define functions

class InvalidConfig(Exception):
	pass 

wildcard_constraints:
	threshold=r"\d+\.\d+",
	separator=r".{0}|_",
	other_flags=r".{0}|[^0-9]+"  # match empty strings or more flags

FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE = "threshold{threshold}{separator}{other_flags}"
DEFAULT_THRESHOLD = .02

MAX_MEM_MB = 250 * 1000  # 250GB

def qc_plot_outputs():
	output_files = []
	for _, row in BIOSAMPLES_CONFIG.iterrows():
		biosample = row["biosample"]
		threshold = determine_threshold(biosample)
		file = os.path.join(RESULTS_DIR, biosample, "Metrics", f"QCSummary_{determine_filtered_prediction_file_format(threshold, config)}.tsv")
		output_files.append(file)
	return output_files

def _get_run_predictions_hic_params(wildcards):
	hic_file = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_file"]
	hic_type = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_type"]
	hic_resolution = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_resolution"]
	if hic_file:
		return f"--hic_file {hic_file} --hic_type {hic_type} --hic_resolution {hic_resolution}"
	else:
		return "--score_column powerlaw.Score"

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressesed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)

def make_paths_absolute(obj, base_path):
	"""
	Use absolute paths to be compatible with github submodules
	Recursively go through the dictionary and convert relative paths to absolute paths.
	"""
	if isinstance(obj, dict):
		for key, value in obj.items():
			obj[key] = make_paths_absolute(value, base_path)
	elif isinstance(obj, str):
		# We assume all strings are paths. If converting the string
		# to an absolute path results in a valid file, then the str was a path
		new_file = os.path.join(base_path, obj)
		if os.path.exists(new_file):
			return new_file
	return obj

def determine_threshold(biosample):
	# config takes priority
	config_threshold = config["params_filter_predictions"]["threshold"]
	if config_threshold:
		return config_threshold
	biosample_row = BIOSAMPLES_CONFIG[BIOSAMPLES_CONFIG["biosample"] == biosample].iloc[0]
	hic_type = biosample_row["HiC_type"]
	if hic_type == None:
		hic_type = "powerlaw"
	elif hic_type == "avg":
		hic_type = "avg_hic"
	elif hic_type == "hic":
		hic_type = "intact_hic"
	matching_row = ABC_THRESHOLDS[
        (ABC_THRESHOLDS["accessibility"] == biosample_row["default_accessibility_feature"])
        & (ABC_THRESHOLDS["has_h3k27ac"] == bool(biosample_row["H3K27ac"]))
        & (ABC_THRESHOLDS["hic_type"] == hic_type)
    ]
	if len(matching_row) == 0:
		print(f"Threshold not found for biosample: {biosample}. Using default threshold of {DEFAULT_THRESHOLD}")
		threshold = DEFAULT_THRESHOLD
	else:
		threshold = matching_row.iloc[0]["threshold"]
	return threshold

def determine_filtered_prediction_file_format(threshold, config):
	include_self_promoter = config['params_filter_predictions']['include_self_promoter']
	only_expressed_genes = config['params_filter_predictions']['only_expressed_genes']
	if include_self_promoter or only_expressed_genes:
		separator = '_'
		other_flags = []
		if include_self_promoter:
			other_flags.append('self_promoter')
		if only_expressed_genes:
			other_flags.append('only_expr_genes')
		other_flags = "__".join(other_flags)
	else:
		separator = ''
		other_flags = ''
	return FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE.format(threshold=threshold, separator=separator, other_flags=other_flags)

def enable_retry(func, func_args={}, max_attempts=3, delay=0.5):
	"""
	To prevent EmptyDataError race condition when using SLURM ro launch jobs as processes
	Assuming the EmptyDataError is caused by a file caching or synchronization lag
	Retry with delay

	@Param
	func:  Function to retry
	func_args:  Dictionary of kwargs for function
	max_attempts:  Maximum number of attempts allowable before raising error
	delay: minimum delay before retry
	"""
	for attempt in range(max_attempts):
		try:
			return func(**func_args)
		except Exception as e:
			if attempt == max_attempts - 1:
				raise
			sleep_time = delay + random.uniform(0, 0.5)
			time.sleep(sleep_time)
	return None

def _validate_accessibility_feature(row: pd.Series):
	if row["DHS"] and row["ATAC"]:
		raise InvalidConfig("Can only specify one of DHS or ATAC for accessibility")
	if not (row["DHS"] or row["ATAC"]):
		raise InvalidConfig("Must provide either DHS or ATAC accessibility file")

def _validate_hic_info(row: pd.Series):
	if row["HiC_file"]:
		if not (row["HiC_type"] and row["HiC_resolution"]):
			raise InvalidConfig("Must provide HiC type and resolution with file")
		# if row["HiC_resolution"] != 5000:
		# 	raise InvalidConfig("Only 5kb resolution supported at the moment")

def _validate_biosamples_config(biosamples_config):
	"""
	Throw exception if a row needs to be fixed
	"""
	for _, row in biosamples_config.iterrows():
		_validate_hic_info(row)
		_validate_accessibility_feature(row)

def load_biosamples_config(config):
	biosamples_config = enable_retry(
		pd.read_csv, 
		func_args={'filepath_or_buffer': config["biosamplesTable"], 'sep': "\t"}
	).replace([np.nan], [None]).set_index("biosample", drop=False)
	biosamples_config["HiC_resolution"] = biosamples_config["HiC_resolution"].replace([None], [0]).astype(int)
	_validate_biosamples_config(biosamples_config)
	_configure_tss_and_gene_files(biosamples_config)
	return biosamples_config

def load_abc_thresholds(config):
	file = config["ref"]["abc_thresholds"]
	return pd.read_csv(file, sep='\t')

def get_accessibility_files(wildcards):
	# Inputs have been validated so only DHS or ATAC is provided
	files = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "DHS"] or BIOSAMPLES_CONFIG.loc[wildcards.biosample, "ATAC"]
	return files.split(",")

def get_peak_files(wildcards):
	files = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "Peak"]
	return files.split(",")

def _configure_tss_and_gene_files(biosamples_config):
	## get TSS and genefile names for each biosample 
	TSS_files = []
	gene_files = []
	for sample in biosamples_config['biosample']:
		tss_file = config['ref']['genome_tss']
		gene_file = config['ref']['genes']
		if biosamples_config.loc[sample, "alt_TSS"]:
			tss_file = biosamples_config.loc[sample, 'alt_TSS']
		if biosamples_config.loc[sample, "alt_genes"]:
			gene_file = biosamples_config.loc[sample, 'alt_genes']
		TSS_files.append(tss_file)
		gene_files.append(gene_file)
					
	biosamples_config["TSS"] = TSS_files
	biosamples_config["genes"] = gene_files

# Making paths absolute is important so that ABC can be 
# used as a submodule for ENCODE-rE2G
# ABC_DIR_PATH = os.path.abspath(config["ABC_DIR_PATH"])
# config = make_paths_absolute(config, ABC_DIR_PATH)
# print(config)

# print(SNAKEFILE_DIR)

RESULTS_DIR = config['results_dir']
BIOSAMPLES_CONFIG = load_biosamples_config(config)
# add estimated powerlaw hic_gamma and hic_scale
existed_cols=BIOSAMPLES_CONFIG.columns.tolist()
if 'hic_gamma' not in existed_cols:
    BIOSAMPLES_CONFIG['hic_gamma']=config['params_predict']['hic_gamma']
if 'hic_scale' not in existed_cols:
    BIOSAMPLES_CONFIG['hic_scale']=config['params_predict']['hic_scale']
#SCRIPTS_DIR = os.path.join(ABC_DIR_PATH, "workflow/scripts")
SCRIPTS_DIR = os.path.join(SNAKEFILE_DIR,"scripts")
ABC_THRESHOLDS = load_abc_thresholds(config)

# rules
# master rule
rule all:
	input:
		allPutative = expand(
			os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"), biosample=BIOSAMPLES_CONFIG["biosample"]
		),
		qcPlots = qc_plot_outputs()

## call macs2 -- if multiple accessibility inputs for one biosample, will aggregate into one output, config['params_macs']['run_macs']
rule call_macs_peaks: 
	input:
		accessibility = get_accessibility_files,
	params:
		pval = config['params_macs']['pval'],
		genome_size = config['params_macs']['genome_size'],
	output: 
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak")
	resources:
		mem_mb=determine_mem_mb
	shell: 
		"""
if [[ "{input.accessibility}" == *tagAlign* ]]; then
	FORMAT="BED"
else
	FORMAT="AUTO"
fi

macs2 callpeak \
-f $FORMAT \
-g {params.genome_size} \
-p {params.pval} \
-n macs2 \
--shift -75 \
--extsize 150 \
--nomodel \
--keep-dup all \
--call-summits \
--outdir {RESULTS_DIR}/{wildcards.biosample}/Peaks \
-t {input.accessibility} 
"""

rule generate_chrom_sizes_bed_file:
	input:
		chrom_sizes = config['ref']['chrom_sizes']
	output:
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed')
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
		awk 'BEGIN {{OFS="\t"}} {{if (NF > 0) print $1,"0",$2 ; else print $0}}' {input.chrom_sizes} > {output.chrom_sizes_bed}
		"""

## sort narrowPeaks
rule sort_narrowpeaks:
	input:
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak") if config['params_macs']['run_macs'] else get_peak_files,
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed')
	params:
		chrom_sizes = config['ref']['chrom_sizes']
	output:
		narrowPeakSorted = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted")
	resources:
		mem_mb=determine_mem_mb
	run:
		if config['params_macs']['run_macs'] or not config['params_macs']['has_header']: ## intersect to remove alternate chromosomes
			shell(f"bedtools intersect -u -a {input.narrowPeak} -b {input.chrom_sizes_bed} | bedtools sort -faidx {params.chrom_sizes} -i stdin > {output.narrowPeakSorted}")
		else:
			shell(f"sed '1d' {input.narrowPeak} | bedtools intersect -u -a stdin -b {input.chrom_sizes_bed} | bedtools sort -faidx {params.chrom_sizes} -i stdin > {output.narrowPeakSorted}")


rule make_candidate_regions:
	input:
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted"),
		accessibility = get_accessibility_files,
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed'),
	params:
		TSS = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'TSS'],
		chrom_sizes = config['ref']['chrom_sizes'],
		regions_blocklist = config['ref']['regions_blocklist'],
		peakExtendFromSummit = config['params_candidate']['peakExtendFromSummit'],
		nStrongestPeak = config['params_candidate']['nStrongestPeaks'],
		output_dir = os.path.join(RESULTS_DIR, "{biosample}", "Peaks"),
		scripts_dir = SCRIPTS_DIR,
		ignoreSummits="--ignoreSummits" if config['params_candidate']['ignoreSummits'] else ""
	output: 
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed")
	resources:
		mem_mb=determine_mem_mb
	shell: 
		"""
python {params.scripts_dir}/makeCandidateRegions.py \
--narrowPeak {input.narrowPeak} \
--accessibility {input.accessibility} \
--outDir {params.output_dir} \
--chrom_sizes {params.chrom_sizes} \
--chrom_sizes_bed {input.chrom_sizes_bed} \
--regions_blocklist {params.regions_blocklist} \
--regions_includelist {params.TSS} \
--peakExtendFromSummit {params.peakExtendFromSummit} \
--nStrongestPeak {params.nStrongestPeak} {params.ignoreSummits}
"""

rule create_neighborhoods:
	input:
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"),
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed')
	params:
		DHS = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "DHS"] or '',
		ATAC = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "ATAC"] or '',
		H3K27ac = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K27ac"] or '',
		default = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'default_accessibility_feature'],
		expression_table = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "Expression"] or '',
		genes = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'genes'],
		ubiquitous_genes = config['ref']['ubiquitous_genes'],
		chrom_sizes = config['ref']['chrom_sizes'],
		qnorm = f"--qnorm {config['ref']['qnorm']}" if config['params_neighborhoods']['use_qnorm'] else "",
		scripts_dir = SCRIPTS_DIR
	output: 
		enhList = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "EnhancerList.txt"),
		geneList = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "GeneList.txt"),
		neighborhoodDirectory = directory(os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods")),
		processed_genes_file = os.path.join(RESULTS_DIR, "{biosample}", "processed_genes_file.bed"),
	resources:
		mem_mb=32*1000
	shell:
		"""
# get sorted & unique gene list
# intersect first to remove alternate chromosomes
bedtools intersect -u -a {params.genes} -b {input.chrom_sizes_bed} | \
bedtools sort -faidx {params.chrom_sizes} -i stdin | \
uniq > {output.processed_genes_file}
						
python {params.scripts_dir}/run.neighborhoods.py \
--candidate_enhancer_regions {input.candidateRegions} \
--default_accessibility_feature {params.default} \
--chrom_sizes {params.chrom_sizes} \
--chrom_sizes_bed {input.chrom_sizes_bed} \
--outdir {output.neighborhoodDirectory} \
--genes {output.processed_genes_file} \
--expression_table {params.expression_table} \
--ubiquitously_expressed_genes {params.ubiquitous_genes} \
{params.qnorm} --DHS {params.DHS} --ATAC {params.ATAC} --H3K27ac {params.H3K27ac}
"""

rule create_predictions:
	input:
		enhancers = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "EnhancerList.txt"),
		genes = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "GeneList.txt"),
	params:
		cellType = lambda wildcards: wildcards.biosample, 
		output_dir = lambda wildcards: os.path.join(RESULTS_DIR, wildcards.biosample, "Predictions"),
		score_column = config['params_filter_predictions']['score_column'],
		hic_params = _get_run_predictions_hic_params,
		chrom_sizes = config['ref']['chrom_sizes'],
		flags = config['params_predict']['flags'],
		# gamma = config['params_predict']['hic_gamma'],
		# scale = config['params_predict']['hic_scale'],
		gamma = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "hic_gamma"],
		scale = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "hic_scale"],
		hic_pseudocount_distance = config['params_predict']['hic_pseudocount_distance'],
		accessibility_feature = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'default_accessibility_feature'],
		scripts_dir = SCRIPTS_DIR,
	output: 
		allPutative = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		allPutativeNonExpressed = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz"),
	resources:
		mem_mb=partial(determine_mem_mb, min_gb=20)  # Use 100GB if using average HiC
	shell:
		"""
python {params.scripts_dir}/predict.py \
--enhancers {input.enhancers} \
--outdir {params.output_dir} \
--score_column {params.score_column} \
--chrom_sizes {params.chrom_sizes} \
--accessibility_feature {params.accessibility_feature} \
--cellType {params.cellType} \
--genes {input.genes} \
--hic_gamma {params.gamma} \
--hic_scale {params.scale} \
--hic_pseudocount_distance {params.hic_pseudocount_distance} \
{params.flags} {params.hic_params}
"""

rule filter_predictions:
	input: 
		allPutative = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		allPutativeNonExpressed = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz"),
	params:
		score_column = config['params_filter_predictions']['score_column'],
		threshold = lambda wildcards: determine_threshold(wildcards.biosample),
		include_self_promoter = config['params_filter_predictions']['include_self_promoter'],
		only_expressed_genes = config['params_filter_predictions']['only_expressed_genes'],
	output:
		enhPredictionsFull = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictionsFull_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
		enhPredictionsFullBedpe = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictionsFull_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.bedpe.gz"),
		enhPredictionsSlim = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictions_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
		genePredictionsStats = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"GenePredictionStats_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv")
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
python {SCRIPTS_DIR}/filter_predictions.py \
--output_tsv_file {output.enhPredictionsFull} \
--output_slim_tsv_file {output.enhPredictionsSlim} \
--output_bed_file {output.enhPredictionsFullBedpe} \
--output_gene_stats_file {output.genePredictionsStats} \
--pred_file {input.allPutative} \
--pred_nonexpressed_file {input.allPutativeNonExpressed} \
--score_column {params.score_column} \
--threshold {params.threshold} \
--include_self_promoter {params.include_self_promoter} \
--only_expressed_genes {params.only_expressed_genes}
"""

rule generate_qc_plot_and_summary:
	input: 
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"),
		neighborhoodDirectory = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods"),
		enhPredictionsFull = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictionsFull_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
		chrom_sizes = config['ref']['chrom_sizes'],
	params:
		output_dir = os.path.join(RESULTS_DIR, "{biosample}", "Metrics"),
		scripts_dir = SCRIPTS_DIR,
		gamma = config['params_predict']['hic_gamma'],
		scale = config['params_predict']['hic_scale'],
	output:
		qc_summary = os.path.join(RESULTS_DIR, "{biosample}", "Metrics", f"QCSummary_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
		qc_plots = os.path.join(RESULTS_DIR, "{biosample}", "Metrics", f"QCPlots_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.pdf")
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
python {params.scripts_dir}/grabMetrics.py \
--outdir {params.output_dir} \
--output_qc_summary {output.qc_summary} \
--output_qc_plots {output.qc_plots} \
--macs_peaks {input.candidateRegions} \
--neighborhood_outdir {input.neighborhoodDirectory} \
--preds_file {input.enhPredictionsFull} \
--chrom_sizes {input.chrom_sizes} \
--hic_gamma {params.gamma} \
--hic_scale {params.scale} 
"""
        
