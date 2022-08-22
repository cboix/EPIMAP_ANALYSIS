#!/usr/bin/env python3
# Author: Benjamin T. James
import h5py
import numpy as np
import pandas as pd
import time, sys

class extractor:
	def __init__(self, pos_mark_list, neg_mark_list, log_file=sys.stderr):
		"""Set shape to Num_locations x n_marks
		   so when index is picked, data is pre-formatted"""
		self.log = log_file
		self.pos_file_list = [h5py.File(f, 'r') for f in pos_mark_list]
		self.neg_file_list = [h5py.File(f, 'r') for f in neg_mark_list]
		for p, n in zip(self.pos_file_list, self.neg_file_list):
			print("Shape: pos=", p['matrix'].shape, file=self.log)
			print("Shape: neg=", n['matrix'].shape, file=self.log)

	def featurize(self, pos_slice, neg_slice, extra_features=None):
		"""Collect indices which match those from the positive selection, in the same number.
		   This confirms that the same DHS indices are sent to both positive and negative
		   in the same number and order."""
		neg_idx = []
		for idx, row in pos_slice.iterrows():
			I = neg_slice.loc[neg_slice['mind'] == row['mind']].sample(n=1).index[0]
			neg_idx.append(I)
		### Unique sorts and removes duplicates with a way to revert to original
		pos_uniq, pos_inverse = np.unique(pos_slice.index, return_inverse=True)
		neg_uniq, neg_inverse = np.unique(neg_idx, return_inverse=True)

		f_pos = np.asarray([p['matrix'][pos_uniq, 0] for p in self.pos_file_list]).reshape(-1, len(self.pos_file_list))
		f_neg = np.asarray([n['matrix'][neg_uniq, 0] for n in self.neg_file_list]).reshape(-1, len(self.neg_file_list))
		f_pos[np.isnan(f_pos)] = 0
		f_neg[np.isnan(f_neg)] = 0
		return [f_pos[pos_inverse, :], f_neg[neg_inverse, :], pos_slice]

	def extract(self, all_info):
		"""Input is a list of several bins,
		the middle one will be the testing bin.
		Format of each bin is the output of featurize()
		"""
		labels = []
		train = []
		test = []
		test_info = None
		### This works since zero-indexing helps the truncated
		### division
		median_index = range(len(all_info))[len(all_info)//2]
		for i, info in enumerate(all_info):
			cur = []
			c_labels = []
			for pos in info[0]:
				cur.append(pos)
				c_labels.append(1)
			for neg in info[1]:
				cur.append(neg)
				c_labels.append(-1)
			# for pos, neg in zip(info[0], info[1]):
			# 	cur.append(pos)
			# 	c_labels.append(1)
			# 	cur.append(neg)
			# 	c_labels.append(-1)
			if i == median_index:
				test = info[0] # just the positive
				test_info = info[2]
			train += cur
			labels += c_labels
		return (np.asarray(train),
				np.asarray(labels),
				np.asarray(test),
				test_info)

	def fprint(self, feat_info, additional_info, cor, cutoff, file_desc):
		df = feat_info.assign(cor=cor)
		threshold = cutoff / (1 + cutoff)
		for i, row in df.loc[df['cor'] >= threshold, :].iterrows():
			out_row = [row['chrom'],
					   row['start'],
					   row['end'],
					   row['gene'],
					   row['cor']] + list(additional_info)
			print("\t".join(map(str, out_row)), file=file_desc)
