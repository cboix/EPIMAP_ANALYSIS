#!/usr/bin/env python3
# Author: Benjamin T. James
import sys
import gzip
import pickle
import numpy as np

if __name__ == "__main__":
	if len(sys.argv[1:]) != 1:
		print("Usage: %s < masterlist.txt masterlist.pickle.gz" % sys.argv[0])
		sys.exit(1)
	with gzip.open(sys.argv[1], 'rb') as P_in:
		indices = pickle.load(P_in)
	index = 0
	masterlist = []
	for rline in sys.stdin:
		L = rline.strip().split()
		masterlist.append(L)
	final = np.asarray(masterlist)[indices]
	for i, row in enumerate(final):
		f_row = list(row[:3])
		f_row.append("%d" % (indices[i]))
		print("\t".join(f_row))
