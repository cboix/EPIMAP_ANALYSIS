#!/usr/bin/env python3
# Author: Benjamin T. James
import sys
import os
import gzip

def for_sample_rna(sample_name, location_table, R, W):
	for rline in R:
		L = rline.decode('utf-8').strip().split()
		if L[0] == sample_name:
			gene = L[1]
			value = float(L[2])
			if gene not in location_table:
				continue
			loc = location_table[gene] # left index, right index
			print("%s\t%d\t%d\t%f\t%s\t%s" % (loc[0], loc[1], loc[2], value, gene, loc[3]), file=W)

def for_sample_rna_val(sample_name, R, W):
	for rline in R:
		L = rline.decode('utf-8').strip().split()
		if L[0] == sample_name:
			print("%s\t%s" % (L[1], L[2]), file=W)

def get_location_table(fname, function=lambda x: x):
	table = {}
	with open(fname, 'r') as R:
		for rline in R:
			line = rline.strip().split()
			if len(line) > 3:
				key = function(line[3])
				table[key] = (line[0], int(line[1]), int(line[2]))
	return table

def for_sample_mark(sample_num, mark_loc_table, R, W):
	for rline in R: #fmt: index sample# value
		L = rline.decode('utf-8').strip().split()
		loc_index = int(L[0])
		sample_index = int(L[1])
		value = float(L[2])
		if sample_index == sample_num and loc_index in mark_loc_table:
			loc = mark_loc_table[loc_index]
			print("%s\t%d\t%d\t%f" % (loc[0], loc[1], loc[2], value), file=W)

def for_sample(sample_name, sample_index, rna_gz, gene_loc_table, mark_loc_table, mark_gz_list):
	print("For sample %s" % sample_name)
	os.mkdir(sample_name)
	with open(os.path.join(sample_name, "rna.txt"), 'w') as W:
		with gzip.open(rna_gz, 'r') as R:
			for_sample_rna(sample_name, gene_loc_table, R, W)
	for mark_gz in mark_gz_list:
		with open(os.path.join(sample_name, os.path.basename(mark_gz) + ".mark"), "w") as W:
			with gzip.open(mark_gz, 'r') as R:
				for_sample_mark(sample_index, mark_loc_table, R, W)

def run(sample_file, rna_gz, gene_location_bed, mark_location_bed, mark_gz_list):
	gene_loc_table = get_location_table(gene_location_bed)
	mark_loc_table = get_location_table(mark_location_bed, function=int)
	with open(sample_file, "r") as R1:
		for rline in R1:
			sample = rline.strip().split()
			sample_name = sample[0]
			sample_num = int(sample[1])
			for_sample(sample_name, sample_num, rna_gz, gene_loc_table, mark_loc_table, mark_gz_list)
if __name__ == "__main__":
	if len(sys.argv[1:]) < 5:
		print("Usage: %s sample_file rna_gz gene_location_bed mark_location_bed mark_1_gz mark_2_gz .." % sys.argv[0])
		sys.exit(1)
	run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:])
#	for_sample(sys.argv[1], sys.argv[3], sys.argv[4], sys.argv[5:])
