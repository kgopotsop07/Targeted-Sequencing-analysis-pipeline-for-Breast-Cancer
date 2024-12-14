#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess


class VariantCalling:
	def __init__(self, reference_sequence, reference_variant_path, aligned_marked_dup_files, variant_calling_out, variant_annotation_path):
		self.reference_sequence = reference_sequence
		self.reference_variant_path = reference_variant_path
		self.aligned_marked_dup_files = aligned_marked_dup_files
		self.variant_calling_out = variant_calling_out
		self.variant_annotation_path = variant_annotation_path
	def somatic_variant_calling(self):
		path, reference = os.path.split(self.reference_sequence)
		get_resource_reference_base = reference [: reference.rfind("_")]
		get_full_variant_ref_path = os.path.join(self.reference_variant_path, get_resource_reference_base)
		somatic_variants = f'''gatk Mutect2 --java-options '-Xmx6G -Dsamjdk.compression_level=5' --native-pair-hmm-threads 6\
		-R {self.reference_sequence} -I {self.aligned_marked_dup_files}.bam -O {self.variant_calling_out}.somatic.vcf.gz\
		-pon {get_full_variant_ref_path}_somatic1000g_pon.b37.vcf --germline-resource {get_full_variant_ref_path}_af.only.gnomad.vcf\
		 ;gatk FilterMutectCalls --java-options '-Xmx6G -Dsamjdk.compression_level=5' --native-pair-hmm-threads 6\
		-R {self.reference_sequence} -V {self.variant_calling_out}.somatic.vcf.gz\
		-O {self.variant_calling_out}.filtered.somatic.vcf.gz --max-events-in-region 4 --filtering-stats {self.variant_calling_out}.somatic.vcf.gz.stats
		'''
		return somatic_variants

	def germline_variant_calling(self):
		germline_variants = f'''gatk HaplotypeCaller --java-options '-Xmx6G -Dsamjdk.compression_level=5' --native-pair-hmm-threads 6\
		-R {self.reference_sequence} -I {self.aligned_marked_dup_files}.bam -O {self.variant_calling_out}.germline.vcf.gz;\
		gatk VariantFiltration --java-options '-Xmx6G -Dsamjdk.compression_level=5' -R {self.reference_sequence} -V {self.variant_calling_out}.germline.vcf.gz\
		-O {self.variant_calling_out}.filtered.germline.vcf.gz --filter-expression 'QD > 2.00' --filter-name "PASS" --filter-expression 'MQ > 40.00'\
		--filter-name "PASS" --filter-expression 'SOR < 3.000' --filter-name "PASS" --filter-expression 'FS < 60.000' --filter-name "PASS"
		'''
		return germline_variants
class VariantAnnotation:
	def __init__(self, variant_calling_path, variant_annotation_path, annotation_references, reference_sequence_file):
		self.variant_calling_path = variant_calling_path
		self.variant_annotation_path = variant_annotation_path
		self.annotation_references = annotation_references
		self.reference_sequence_file = reference_sequence_file
	@staticmethod
	def _get_reference_base(reference_file):
		path, reference_file_name = os.path.split(reference_file)
		reference_file_suffix = reference_file_name[: reference_file_name.rfind("_")]
		return reference_file_suffix
	def get_somatic_annotations(self):
		get_reference_base = VariantAnnotation._get_reference_base(self.reference_sequence_file)
		full_annotation_reference_path = os.path.join(self.annotation_references,get_reference_base)
		annotated_somatic_variants = f'''table_annovar.pl {self.variant_calling_path}.filtered.somatic.vcf.gz {full_annotation_reference_path} \
		-buildver {get_reference_base} --protocol refGene,cosmic,clinvar_20240611,avsnp147 --operation g,f,f,f \
		--outfile {self.variant_annotation_path}-somatic --remove --nastring . --vcfinput

		'''
		return annotated_somatic_variants
	def get_germline_annotations(self):
		get_reference_base = VariantAnnotation._get_reference_base(self.reference_sequence_file)
		full_annotation_reference_path = os.path.join(self.annotation_references,get_reference_base)
		annotated_germline_variants = f'''table_annovar.pl {self.variant_calling_path}.filtered.germline.vcf.gz {full_annotation_reference_path} \
		-buildver {get_reference_base} --protocol refGene,avsnp147,clinvar_20240611 --operation g,f,f \
		--outfile {self.variant_annotation_path}-germline  --nastring . --vcfinput

		'''
		return annotated_germline_variants


def get_variants(analysis_type, marked_dup_files, reference_file, variant_called_file_base, variant_annotation_file_base, annotation_reference_base_folder, variant_path_reference):
	init_variant_detection = VariantCalling(reference_file, variant_path_reference, marked_dup_files, variant_called_file_base, variant_annotation_file_base)
	init_variant_annotation = VariantAnnotation(variant_called_file_base, variant_annotation_file_base, annotation_reference_base_folder, reference_file)
	if analysis_type == "somatic":
		print (f"=>... accessing somatic variant detected in {marked_dup_files}")
		somatic_variant_command = init_variant_detection.somatic_variant_calling()
		subprocess.run(somatic_variant_command, capture_output = True, shell = True)
		print (f"=>... performing somatic variant annotation in {marked_dup_files}")
		somatic_variant_annotation = init_somatic_annotation.get_somatic_annotations()
		subprocess.run (somatic_variant_annotation, capture_output = True, shell = True)

	elif analysis_type == "germline":
		print (f"=>... accessing germline variant detected in {marked_dup_files}")
		germline_variant_command = init_variant_detection.germline_variant_calling()
		subprocess.run (germline_variant_command, capture_output = True, shell = True)
		print (f"=>... performing germline variant annotation in {marked_dup_files}")
		germline_variant_annotation = init_variant_annotation.get_germline_annotations()
		subprocess.run (germline_variant_annotation, capture_output = True, shell = True)

