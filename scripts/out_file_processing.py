#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess

class PrepSeqAlign:
	def __init__(self, aligned_full_path):
		self.aligned_full_path = aligned_full_path
	def align_files_to_ref(self):
		if self.aligned_full_path.endswith(".sam"):
			align_base = self.aligned_full_path[: self.aligned_full_path.rfind(".sam")]
			full_path_base = os.path.join(self.aligned_full_path, align_base)
			samtools_command = ''' {1} view -@ 6 -S -b {0}.sam > {0}.bam;
			{1} sort -@ 6 -n {0}.bam -o {0}.sort.base.bam;
			{1} fixmate -@ 6 -m {0}.sort.base.bam {0}.fixmate.bam;
			{1} sort -@ 6 {0}.fixmate.bam -o {0}.sorted_file.bam;
			{2} AddOrReplaceReadGroups I={0}.sorted_file.bam O={0}.sorted_reads_with_RG.bam SORT_ORDER=coordinate RGID=4 RGLB=lib2 RGPL=illumina RGSM=30 CREATE_INDEX=True RGPU=unit1 RGPI=300;
			{2} MarkDuplicates I={0}.sorted_reads_with_RG.bam O={0}.untrained_sorted_markedup.bam M={0}.untr_metrics.txt
			'''.format(align_base, "samtools","picard-tools")
			return samtools_command
	@staticmethod
	def get_targeted_regions (targeted_sequencing_file, output_marked_files):
		if targeted_sequencing_file:
			sam_targeted_panel = f"samtools view -b -L {targeted_sequencing_file} {output_marked_files} > {output_marked_files}.targeted.trained_sorted_markedup"
			return sam_targeted_panel

class PrepVarFiles:
	def __init__(self, input_reference_file, gatk_folder):
		self.input_reference_file = input_reference_file
		self.gatk_folder = gatk_folder
	@staticmethod
	def get_vcf_base(vcf_input_path):
		vcf_file_prefix = vcf_input_path[ :vcf_input_path.rfind(".vcf")]
		return vcf_file_prefix
	def get_gatk_path(self):
		full_reference_path_split = os.path.split(self.input_reference_file)
		reference_type = full_reference_path_split[1][ :full_reference_path_split[1].rfind("_")]
		gatk_reference_path = os.path.join(self.gatk_folder,reference_type)
		return gatk_reference_path

def alignment_handling(sam_outfile, input_ref, gatk_res_folder, file_path, reference_bed_file = None):
	print ("=> processing sequence aligned and mapped {.SAM} files")
	init_run_alignment = PrepSeqAlign(sam_outfile)
	run_samtools = init_run_alignment.align_files_to_ref()
	subprocess.run(run_samtools, capture_output = True, shell = True)
	if f"{sam_outfile}".endswith(".sam"):
		os.remove(f"{sam_outfile}.sam")
	else:
		pass
	if reference_bed_file is not None:
		run_targeted_Sequence = PrepSeqAlign.get_targeted_regions(reference_bed_file, sam_outfile) 
		subprocess.run (run_targeted_Sequence, capture_output = True, shell = True)
	else:
		pass
	print ("=> preparing gatk resources files")
	init_run_gatk_preps = PrepVarFiles(input_ref, gatk_res_folder)
	check_gatk_indexes_base = init_run_gatk_preps.get_gatk_path()
	for vcf_file in os.listdir(check_gatk_indexes_base):
		vcf_file_prefix = PrepVarFiles.get_vcf_base(vcf_file)
		vcf_path_and_prefix = os.path.join(check_gatk_indexes_base,vcf_file_prefix)
		if os.path.isfile(vcf_path_and_prefix + ".vcf.idx"):
			pass
		else:
			run_gatk_indexes = f"gatk IndexFeatureFile -I {vcf_path_and_prefix}.vcf"
			subprocess.run (run_gatk_indexes, capture_output = True, shell = True)

	print ("=> preparing the reference genome for variant calling")
	reference_seq_base = input_ref[ :input_ref.rfind(".fasta")]
	if not os.path.isfile(reference_seq_base + ".dict") and not os.path.isfile(input_ref + "fai"):
		prepare_required_reference_files = f"gatk CreateSequenceDictionary -R {input_ref}; samtools faidx {input_ref}"
		subprocess.run(f"{prepare_required_reference_files}", capture_output = True, shell = True)
	else:
		print ("all necessary files are present")

if __name__=="__main__":
	alignment_handling(sam_output_file, ref_input_file, gat_resource_folder, in_file_path)
