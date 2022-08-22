#!/usr/bin/env python3
# Author: Benjamin T. James
import h5py
import numpy as np
import sys
import pandas as pd

class chromhmm_state:
	def __init__(self, fname, log_file=sys.stderr):
		self.fname = fname
		self.fd = None
		self.labels = []
		self.table = {}
		self.log = log_file

	def __enter__(self, encoding='ascii'):
		"""Open the DHS ChromHMM file and
		build associated data structures"""
		self.fd = h5py.File(self.fname, 'r')
		print("[DEBUG] Entering ChromHMM reading from %s" % self.fname, file=self.log)
		self.labels = [lab.decode(encoding) for lab in self.fd['sample_list'][:]]
		print("Labels:", self.labels[:min(0, len(self.labels))], "shape=", self.fd['sample_list'].shape, file=self.log)
		enh_state_list = [lab.decode(encoding) for lab in self.fd["chromhmm_state_list"][:]]

		enh_flag_list = self.fd["chromhmm_flag_list"][:]
		self.dhs_names = np.asarray([lab.decode(encoding) for lab in self.fd['location_list'][:]])
		self.table = {k: v for k, v in zip(enh_state_list, enh_flag_list)}
		print("Enh lookup table:", self.table, file=self.log)
		return self

	def sample_names_to_indices(self, which_labels):
		lab_index_list = []
		for label in which_labels:
			if label not in self.labels:
				raise ValueError("%s is not a valid sample name" % label)
			else:
				lab_index_list.append(self.labels.index(label))
		return lab_index_list

	def get_dhs_names(self, which_labels, which_states):
		"""Get items in selected sample(s), which are columns,
		such that the rows when masked with the enhancer flags
		are nonzero. Return the DHS names which match those indices.
		"""
		if not self.fd:
			raise RuntimeError("Must use __enter__")
		if not which_labels or type(which_labels) != list:
			raise ValueError("which_labels is not a list")
		if not which_states or type(which_states) != list:
			raise ValueError("which_states is not a list")
		lab_index_list = []
		for label in which_labels:
			if label not in self.labels:
				raise ValueError("%s is not a valid sample name" % label)
			else:
				lab_index_list.append(self.labels.index(label))
		print("[DEBUG] Label indices:", which_labels, "->", lab_index_list, file=self.log)
		flag = 0
		for state in which_states:
			if state not in self.table:
				raise ValueError("%s is not a valid state" % state)
			flag |= self.table[state]
		print("[DEBUG] Enhancer flags:", which_states, "->", flag, file=self.log)
		dhs_indices = np.where(self.fd['data'][:,lab_index_list] & flag)[0]
		print("DHS indices:", dhs_indices, file=self.log)
		dhs_names = self.dhs_names[dhs_indices]
		return dhs_names

	def __exit__(self, e_type, e_val, e_tb):
		print("[DEBUG] Exiting ChromHMM file", self.fname, file=self.log)
		self.fd.close()
		self.fd = None

if __name__ == "__main__":
	if len(sys.argv[1:]) != 4:
		print("usage: %s chunks.bed in.h5 sample state")
		sys.exit(1)
	in_bed = pd.read_csv(sys.argv[1], sep="\t", header=None)
	with chromhmm_state(sys.argv[2]) as cs:
		dhs_names = cs.get_dhs_names(which_labels=[sys.argv[3]], which_states=[sys.argv[4]])
	in_bed.loc[in_bed[3].isin(dhs_names)].to_csv(sys.stdout, sep="\t", index=False, header=False)
