#!/usr/bin/env python3
# Author: Benjamin T. James
import sys, os
import requests
import gzip
import pandas as pd
import numpy as np
import pybedtools
import h5py
import argparse as ap

def get_site(sample_name):
	site = "https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg19/CALLS/%s_18_CALLS_segments.bed.gz" % (sample_name.upper())
	r = requests.get(site)
	strng = gzip.decompress(r.content).decode('utf-8')
	return pybedtools.BedTool(strng, from_string=True)

# {0, 1, 2, 3} are enhancer loc. and type
# {4, 5, 6, 7} are masterlist loc and index
def get_matrix(masterlist_bed, sample_list):
	m_df = pd.read_csv(masterlist_bed, sep="\t", header=None)
	# Use int32, max=2^32 > 2^18, which is the largest possible flag due to 18 state ChromHMM
	M = np.zeros((len(m_df), len(sample_list)), dtype='int32')
	# Have a lookup table of ChromHMM states
	# which holds flags which can be combined
	state_tbl = {}
	name_tbl = {n: i for i, n in enumerate(m_df[3].values)}
	for sample_idx, sample in enumerate(sample_list):
		print(sample_idx, sample)
		bed = get_site(sample)

		# wa and wb print both bed files so the data in the last columns are extracted
		res = bed.intersect(masterlist_bed, wa=True, wb=True)
		df = res.to_dataframe()
		os.unlink(bed.fn)
		os.unlink(res.fn)
		# pybedtools does not delete temp files, so delete after use
		# This is also why 'masterlist_bed' is re-read every loop,
		# so it is not accidentally deleted.

		state = df.values[:,3] # col[0-3]=state col[4-7]=masterlist
		name_list = df.values[:,7]
		# Add states if do not exist in master state table
		for s in set(state):
			if s not in state_tbl:
				# Get an information-preserving flag
				# e.g. binary power
				# E12 -> (1 << 12) = 2^12
				val = 1 << int(s[1:])
				state_tbl[s] = val
		for s, name in zip(state, name_list):
			enh_idx = state_tbl[s]
			name_idx = name_tbl[name]
			M[name_idx, sample_idx] |= enh_idx
		print(M[:100,sample_idx])
	return M, state_tbl, m_df[3].values

if __name__ == "__main__":
	if len(sys.argv[1:]) != 3:
		print("Usage: %s masterlist_indices_with_chunk_names sample_list out.h5" % sys.argv[0])
		print("\twhere sample_list is a line separated file of sample names")
		sys.exit(1)
	sample_list = pd.read_csv(sys.argv[2], sep="\t", header=None).values[:,0]
	M, state_tbl, name_list = get_matrix(sys.argv[1], sample_list)

	state_list = list(state_tbl.keys())
	# HDF5 format must have fixed-length strings
	enc_state_list = [np.string_(s) for s in state_list]
	enc_name_list = [np.string_(s) for s in name_list]
	state_flags = [state_tbl[s] for s in state_list]

	h = h5py.File(sys.argv[3], 'w')
	h.create_dataset('data', data=M, compression="gzip", chunks=True)

	enc_sample_list = [np.string_(s) for s in sample_list]
	h.create_dataset('location_list', data=enc_name_list, compression="gzip", chunks=True)
	h.create_dataset('sample_list', data=enc_sample_list, compression="gzip", chunks=True)
	# must print out chromhmm data so mat. can be decoded
	h.create_dataset('chromhmm_state_list', data=enc_state_list)
	h.create_dataset('chromhmm_flag_list', data=state_flags)
	h.close()
