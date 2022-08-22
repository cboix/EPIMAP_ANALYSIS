#!/usr/bin/env python3
# Author: Benjamin T. James
import sys
import os
import gzip
from sample_grep import for_sample_rna_val

def for_sample(sample_name, rna_gz):
	print("For sample %s" % sample_name)
	if not os.path.isdir(sample_name):
		os.mkdir(sample_name)
	with open(os.path.join(sample_name, "rna.txt"), 'w') as W:
		with gzip.open(rna_gz, 'r') as R:
			for_sample_rna_val(sample_name, R, W)
	return 0

def run(sample_file, rna_gz):
#	gene_loc_table = get_location_table(gene_location_bed)
	with open(sample_file, "r") as R1:
		for rline in R1:
			sample = rline.strip().split()
			sample_name = sample[0]
			sample_num = int(sample[1])
			for_sample(sample_name, rna_gz)
	return 0

if __name__ == "__main__":
	if len(sys.argv[1:]) != 2:
		print("Usage: %s sample_file rna_gz" % sys.argv[0])
		sys.exit(1)
	run(sys.argv[1], sys.argv[2])
#	for_sample(sys.argv[1], sys.argv[3], sys.argv[4], sys.argv[5:])
