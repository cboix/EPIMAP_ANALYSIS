#!/usr/bin/env python3
# Author: Benjamin T. James
import sys
from get_tss import get_master_chrom_list
import numpy as np

class Bed:
	def __init__(self, chrom_size_table, chrom_order, interval):
		self.chrom_order = chrom_order
		self.interval = interval

		# creates vectors corresponding to length of each chrom
		self.n_bins = [int(np.ceil(chrom_size_table[k] / interval)) for k in chrom_order]
		self.values = [np.zeros((v,)) for v in self.n_bins]
		self.count = [np.zeros((v,)) for v in self.n_bins]

	def add(self, chrom, left, right, value=1):
		l_bin = left // self.interval
		r_bin = (right-1) // self.interval # end not inclusive
		idx = self.chrom_order.index(chrom)
		for n in range(l_bin, r_bin+1):
			self.values[idx][n] += value
			self.count[idx][n] += 1
		return 0

	def get(self):
		out = {}
		for i, chrom in enumerate(self.chrom_order):
			count = self.count[i]
			count[count == 0] = 1
			out[chrom] = self.values[i] / count
		return out


def get_bed(chrom_size, interval, R):
	mc_list = get_master_chrom_list()
	b = Bed(chrom_size, mc_list, interval)
	for r_line in R:
		line = r_line.strip().split()
		chrom = line[0]
		begin = int(line[1])
		end = int(line[2])
		if len(line) > 3:
			val = float(line[3])
		else:
			val = 1
		b.add(chrom, begin, end, val)
	return b.get()


def get_chrom_size_table(fname):
	table = {}
	with open(fname, 'r') as chromR:
		for rline in chromR:
			L = rline.strip().split()
			table[L[0]] = int(L[1])
	return table


def run(chrom_size, interval, R, W):
	bed = get_bed(chrom_size, interval, R)


if __name__ == "__main__":
	if len(sys.argv[1:]) != 3:
		print("Usage: %s bed chrom.sizes block_size" % sys.argv[0])
		sys.exit(1)
	table = get_chrom_size_table(sys.argv[2])
	with open(sys.argv[1], 'r') as R:
		run(table, int(sys.argv[3]), R, sys.stdout)
