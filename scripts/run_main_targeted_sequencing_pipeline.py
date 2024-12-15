#!/usr/bin/env python3
## -----------------------------------------------------------------------------------
## Author	: Kgopotso Phakwago
## Description	: Full workflow for analysis of NGS targeted sequencing data analysis
## -----------------------------------------------------------------------------------

import os
import sys
import argparse
import subprocess
from data_alignment import execute_alignment_step


def run_targeted_sequence_pipeline():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input_path", help = "please specify the paired-end fastQ files with extension <{file}_1/2.fastq> or <{input}_R1/R2_001.fastq>")
	parser.add_argument("-r", "--reference_sequence_file", help = "specify the human reference sequence.fasta (not fna). Remember the prefix starts with the reference type <hg19 or hg38>")
	parser.add_argument("-p", "--bqsr", help = "perform base quality score recalibration (bqsr)?", choices = ["perform", "ignore"])
	parser.add_argument("-a", "--analysis_type", help = "Perform analysis for <germline variants> or <somatic variants> ?",choices = ["somatic", "germline"] )
	parser.add_argument("-b", "--targeted_regions", help = "please specify the targeted regions used for preparing the illumina panel libraries, Optional", default = None)
	args = parser.parse_args()
	##------------SPECIFY INPUT DB PATH------------------
	WD = os.getcwd()
	databases = os.path.split(WD)[0]
	##-------------------------------------------------
	get_variant_reference_file_path = os.path.join(databases,"reference_files","variant_calling_reference_files")
	get_annotation_reference_file_path = os.path.join(databases,"reference_files","annotationt_reference_files")
	get_reference_prefix = os.path.split(args.reference_sequence_file)[1][: os.path.split(args.reference_sequence_file)[1].rfind("_")]
	##------------------RUN PIPELINE-------------------
	if args.input_path is None or args.reference_sequence_file is None:
		subprocess.run('bash -c "python3 run_main_targeted_sequencing_pipeline.py -h"', capture_output = True, shell = True)
	else:
		execute_alignment_step(args.input_path, args.reference_sequence_file, f"{get_variant_reference_file_path}", args.bqsr, args.analysis_type, f"{get_annotation_reference_file_path}", args.targeted_regions)
		subprocess.run (f"Rscript {WD}/get_combined_maf_table.R {get_annotation_reference_file_path} {get_reference_prefix}", capture_output = True, shell = True)

if __name__=="__main__":
	run_targeted_sequence_pipeline()




