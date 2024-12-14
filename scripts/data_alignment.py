#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from create_directories import create_required_dir
from out_file_processing import alignment_handling
from try_base_quality_recalibration import perform_base_quality_recalibration
from perform_variant_calling_and_annotation import get_variants


class GetFastqBase:
	def __init__(self, suffix,in_file):
		self.__suffix = suffix
		self.__in_file = in_file
	def get_file_base(self):
		out_file_name = self.__in_file[: self.__in_file.rfind(str(self.__suffix))]
		return out_file_name
class ReferenceGnome:
	def __init__(self, GRCh_ref):
		self.GRCh_ref = GRCh_ref
	def select_reference_genome_path(self):
		ref_path, human_ref = os.path.split(self.GRCh_ref)
		reference_prefix = human_ref[: human_ref.rfind("_")]
		return ref_path, human_ref, reference_prefix, self.GRCh_ref
	@staticmethod
	def get_bwa_index(ref_file_name):
		ref_variant = ReferenceGnome(ref_file_name)
		get_ref_file_name = ref_variant.select_reference_genome_path() [3]
		if os.path.isfile(str(get_ref_file_name)+".bwt"):
			print ("Burrows Wheler index files exist in {0}".format(os.path.split(ref_file_name)[0]))
		else:
			subprocess.run(f"bash -c 'bwa index {get_ref_file_name}'", capture_output = True, shell = True)
def execute_alignment_step(file_path, ref, res_bundle, recalibration, analysis, annotation_reference_base_fol, targeted_bed_file = None):
	init_ref_variant = ReferenceGnome(ref)
	ref_folder = init_ref_variant.select_reference_genome_path() [3]
	ReferenceGnome.get_bwa_index(ref)
	variant_path_ref = os.path.join(res_bundle, ref[: ref.rfind("_")])
	for folder_list in os.listdir(file_path):
		if "output_directory" not in folder_list:
			create_required_dir(file_path)
		else:
			pass
	for files in os.listdir(file_path):
		if files.endswith("_1.fastq.gz"):
			init_file_base = GetFastqBase("_1.fastq.gz", files)
			fastq_base = init_file_base.get_file_base()
			sam_out_put = os.path.join(file_path,"output_directory","read_aligned_files",fastq_base)
			variant_output_base = os.path.join(file_path,"output_directory","variant_called_files",fastq_base)
			variant_annotation_base = os.path.join(file_path,"output_directory","variant_annotated_files",fastq_base)
			file_full_path = os.path.join(file_path,fastq_base)
			subprocess.run (f"bash -c 'bwa mem -t 6 -k 30 -M {ref_folder} {file_full_path}_1.fastq.gz {file_full_path}_2.fastq.gz > {sam_out_put}.sam '", capture_output = True, shell = True)
			alignment_handling(f"{sam_out_put}.sam", ref, res_bundle, file_path)
			if recalibration == "perform":
				perform_base_quality_recalibration (ref, res_bundle, f"{sam_out_put}")
				get_variants(analysis, f"{sam_out_put}.trained_sorted_markedup", ref, variant_output_base, variant_annotation_base, annotation_reference_base_fol, res_bundle)
			elif recalibration == "ignore":
				subprocess.run(f"bash -c 'samtools index {sam_out_put}.untrained_sorted_markedup.bam'", capture_output = True, shell = True)
				get_variants(analysis, f"{sam_out_put}.untrained_sorted_markedup", ref, variant_output_base, variant_annotation_base, annotation_reference_base_fol, res_bundle)
		elif files.endswith("_R1_001.fastq"):
			init_file_base = GetFastqBase("_R1_001.fastq.gz", files)
			fastq_base = init_file_base.get_file_base()
			sam_out_put = os.path.join(file_path,"output_directory","read_aligned_files",fastq_base)
			variant_output_base = os.path.join(file_path,"output_directory","variant_called_files",fastq_base)
			variant_annotation_base = os.path.join(file_path,"output_directory","variant_annotated_files",fastq_base)
			file_full_path = os.path.join(file_path,fastq_base)
			subprocess.run (f"bash -c'bwa mem -t 6 -k 30 -M {ref_folder} {file_full_path}_R1_001.fastq.gz {file_full_path}_R2_001.fastq.gz > {sam_out_put}.sam '", capture_output = True, shell = True)
			alignment_handling(f"{sam_out_put}.sam", ref, res_bundle, file_path, targeted_bed_file)
			if recalibration == "perform":
				perform_base_quality_recalibration (ref, res_bundle, f"{sam_out_put}")
				get_variants(analysis, f"{sam_out_put}.trained_sorted_markedup", ref, variant_output_base, variant_annotation_base, annotation_reference_base_fol, res_bundle)
			elif recalibration == "ignore":
				subprocess.run(f"bash -c 'samtools index {sam_out_put}.untrained_sorted_markedup.bam'", capture_output = True, shell = True)
				get_variants(analysis, f"{sam_out_put}.untrained_sorted_markedup", ref, variant_output_base, variant_annotation_base, annotation_reference_base_fol, res_bundle)
		else:
			sys.stderr.write('please provide fastq paired-end files with either sufix {prefix}_1/2.fastq.gz or {prefix}_R1/R2_001.fastq.gz\n')


if __name__=="__main__":
	execute_alignment_step(input_dirs, ref_file, gatk_res_path, base_qual_recal, variants_analysis_type, annotation_references, bed_file)
