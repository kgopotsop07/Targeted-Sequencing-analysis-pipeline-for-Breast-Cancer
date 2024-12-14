#!/usr/bin/env python3

import os
import sys
import re
import Bio
from Bio import SeqIO
args = sys.argv[1]
## ----- adding arguments ------
list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","X","Y"]
def get_input_files(file_name):
	fh = SeqIO.parse(open(file_name, "r"), "fasta")
	fw = open("Homo_sapiens_assembly_GRCh38_Gatk_latest_genomic.ensembl.fna", "w+")
	for record in fh:
		sequences = record.seq
		entries = record.id
		if entries in list:
			Z = ">chr"+str(entries)+"\n"+str(sequences)+"\n"
			fw.write(Z)
get_input_files(args)
