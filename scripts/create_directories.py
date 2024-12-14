#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess

def create_required_dir(input_dir):
	out_base_folder = "output_directory"
	out_dirs = ["read_aligned_files", "variant_annotated_files", "variant_called_files"]
	for out_dir in out_dirs:
		folder = os.path.join(input_dir,out_base_folder,out_dir)
		if os.path.exists(folder):
			pass
		else:
			os.makedirs(folder)
if __name__=="__main__":
	create_required_dir(fastq_dirs)

