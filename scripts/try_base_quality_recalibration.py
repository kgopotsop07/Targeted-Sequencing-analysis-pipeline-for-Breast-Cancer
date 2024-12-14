#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess

class PrepareVariantReference:
	def __init__(self, reference_sequence, variant_resource_folder):
		self.__reference_sequence = reference_sequence
		self.variant_resource_folder = variant_resource_folder
	@staticmethod
	def __get_reference_base(reference_file_path):
		path, reference_file = os.path.split(reference_file_path)
		reference_base = reference_file[: reference_file.rfind("_")]
		return reference_base
	def train_variant_resource_files(self):
		reference_sequence_prefix = PrepareVariantReference.__get_reference_base(self.__reference_sequence)
		variant_resource_path = os.path.join(self.variant_resource_folder, reference_sequence_prefix)
		return variant_resource_path

def perform_base_quality_recalibration(human_reference_file, resource_bundle_variant_path, untr_bam_file_base):

	init_path = PrepareVariantReference(human_reference_file, resource_bundle_variant_path)
	variant_reference_file_path = init_path.train_variant_resource_files()
	list_of_variant_refs = []
	for variant_reference_files in os.listdir(variant_reference_file_path):
		varinat_full_path = os.path.join(variant_reference_file_path, variant_reference_files)
		list_of_variant_refs.append(varinat_full_path)
	condition = lambda x: x.endswith("vcf")
	filtered_list_of_variant_refs = [x for x in list_of_variant_refs if condition(x)]
	gatk_base_qual_recal = '''gatk BaseRecalibrator --java-options '-Xmx6G -Dsamjdk.compression_level=5'\
	-I {0}.untrained_sorted_markedup.bam -O {0}.1_recal_data.table -R {1} --known-sites {2}\
	--known-sites {3} --known-sites {4}; gatk ApplyBQSR --java-options '-Xmx6G -Dsamjdk.compression_level=5'\
	-bqsr {0}.1_recal_data.table -I {0}.untrained_sorted_markedup.bam -O {0}.trained_sorted_markedup.bam
	'''.format(untr_bam_file_base, human_reference_file, filtered_list_of_variant_refs [0], filtered_list_of_variant_refs [1], filtered_list_of_variant_refs[2])
	subprocess.run(gatk_base_qual_recal, capture_output = True, shell = True)
