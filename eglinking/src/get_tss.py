#!/usr/bin/env python3
# Author: Benjamin T. James
import sys
import gzip
import numpy as np

master_chrom_list = sorted(list(map(lambda x: "chr%d" % x, range(1, 23))) + ["chrX"])

def get_master_chrom_list():
	return list(master_chrom_list)

def get_tss(gene_info_gz):
	table = {}
	with gzip.open(gene_info_gz, 'r') as ginfo:
		for line in ginfo:
			L = line.decode('utf8').strip().split()
			if L[0][0] == "#":
				continue
			if len(L) > 9 and L[2] == "transcript" and L[0] in master_chrom_list:
				chrom = L[0]
				left = int(L[3])
				right = int(L[4])
				name = L[9]
				strand = L[6]
				if strand == '+':
					tss = left
				else:
					tss = right
				if name in table:
					tval = table[name]
					if tval[2] == '+' and tss < tval[1]:
						table[name] = (chrom, tss, strand)
					elif tval[2] == '-' and tss > tval[1]:
						table[name] = (chrom, tss, strand)
				else:
					table[name] = (chrom, tss, strand)
	return table

def get_uniq_filter_list(gz_name):
	gset = set()
	with gzip.open(gz_name, 'r') as R:
		for rline in R:
			gname = rline.decode('utf').strip().split()[1]
			if gname[:3] == "ENS":
				gset.add(gname)
	return gset

def run(gene_info_gz, f_out, filter_gz=None):
	table = get_tss(gene_info_gz)
	gset = set()
	if filter_gz:
		gset = get_uniq_filter_list(filter_gz)
	for k, v in table.items():
		gene = k.split('.')[0].replace('"', '').replace(';', '') # remove version number (dot suffix), quotes, and semicolons
		chrom = v[0]
		tss = v[1]
		strand = v[2]
		if not(filter_gz) or gene in gset:
			print("%s\t%s\t%d\t%s" % (gene, chrom, tss, strand), file=f_out)

if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("Usage: %s gene_info.gz [filter.gz]" % sys.argv[0])
		sys.exit(1)
	if len(sys.argv) == 2:
		run(sys.argv[1], sys.stdout)
	else:
		run(sys.argv[1], sys.stdout, filter_gz=sys.argv[2])
