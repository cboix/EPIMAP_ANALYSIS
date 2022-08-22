#!/usr/bin/env python3
# Author: Benjamin T. James
# Modified by Carles Boix
import h5py
import numpy as np
import pandas as pd
import time, sys
from collections import defaultdict
import gc

class feature:
    def __init__(self, pos_mark_list, neg_mark_list, log_file=sys.stderr):
        """Set shape to Num_locations x n_marks
           so when index is picked, data is pre-formatted"""
        self.log = log_file
        print("Initializing feature object")
        print("Positive files:", pos_mark_list)
        print("Negative files:", neg_mark_list)
        self.pos_file_list = [h5py.File(f, 'r') for f in pos_mark_list]
        self.neg_file_list = [h5py.File(f, 'r') for f in neg_mark_list]
        for p, n in zip(self.pos_file_list, self.neg_file_list):
            print("Shape: pos=", p['matrix'].shape, file=self.log)
            print("Shape: neg=", n['matrix'].shape, file=self.log)

    def gen_length_feature(min_length, max_length, exponent=-1):
        # offset is not inside the N feature,
        # so must encode from is positive other
        return lambda p, n: (np.abs(p['offset'])**exponent,
                             np.abs(p['offset'])**exponent)

    def featurize(self, pos_slice, neg_slice, extra_features=[], combine=False):
        """Collect indices which match those from the positive selection, in the same number.
        This confirms that the same DHS indices are sent to both positive and negative
        in the same number and order."""
        print("Obtaining features from hdf5")
        # Matching negative indices (very slow)
        print("Matching negative indices:")
        t1 = time.time()
        pos_slice['chrom_mind'] = pos_slice['chrom'] + "_" + \
            pos_slice['mind'].map(str)
        neg_slice['chrom_mind'] = neg_slice['chrom'] + "_" + \
            neg_slice['mind'].map(str)
        pos_slice_uq = pos_slice[['chrom_mind']].drop_duplicates()
        # Turn into dictionary:
        print("Making positive dictionary")
        pos_dict = defaultdict(list)
        for key, value in zip(pos_slice['chrom_mind'],
                              list(range(pos_slice.shape[0]))):
            pos_dict[key].append(value)
        print("Making negative dictionary")
        neg_dict = defaultdict(list)
        for key, value in zip(neg_slice['chrom_mind'], neg_slice.index):
            neg_dict[key].append(value)
        i = 0
        neg_idx = np.zeros(pos_slice.shape[0], dtype=int)
        for key, pindices in pos_dict.items():
            i += 1
            npos = len(pindices)
            I = np.random.choice(neg_dict[key], size=npos, replace=True)
            neg_idx[pindices] = I
        print(time.time() - t1)
        # Old version
        # neg_idx = []
        # for idx, row in pos_slice.iterrows():
        #     I = neg_slice.loc[neg_slice['mind'] == row['mind']].sample(n=1).index[0]
        #     neg_idx.append(I)
        ### Unique sorts and removes duplicates with a way to revert to original
        print("Get unique positions:")
        pos_uniq, pos_inverse = np.unique(pos_slice.index, return_inverse=True)
        neg_uniq, neg_inverse = np.unique(neg_idx, return_inverse=True)
        # Initialize positive and negative feature matrices:
        print("Initializing positive and negative feature matrices:")
        f_pos = np.zeros((len(pos_uniq), len(self.pos_file_list) + len(extra_features)))
        f_neg = np.zeros((len(neg_uniq), len(self.neg_file_list) + len(extra_features)))
        print("Postive features matrix:", f_pos.shape)
        print("Negative features matrix:", f_neg.shape)
        # Populate from the hdf5 files:
        for i, (pfile, nfile) in enumerate(zip(self.pos_file_list,
                                               self.neg_file_list)):
            print('Loading feature set: ', i)
            pmat = pfile['matrix'][:,0]
            f_pos[:, i] = pmat[pos_uniq]
            del(pmat)
            gc.collect()
            nmat = nfile['matrix'][:,0]
            f_neg[:, i] = nmat[neg_uniq]
            del(nmat)
            gc.collect()
            # f_pos[:, i] = pfile['matrix'][pos_uniq, 0]
            # f_neg[:, i] = nfile['matrix'][neg_uniq, 0]
        # For extra features
        for i, feat in enumerate(extra_features):
            print('Loading extra features: ', i)
            (pdata, ndata) = feat(pos_slice, neg_slice.loc[neg_idx, :])
            f_pos[pos_inverse, len(self.pos_file_list) + i] = pdata
            f_neg[neg_inverse, len(self.neg_file_list) + i] = ndata
        f_pos[np.isnan(f_pos)] = 0
        f_neg[np.isnan(f_neg)] = 0
        print("Getting features took:", time.time() - t1)
        if combine:
            normal_L = len(self.pos_file_list)
            for j in range(len(extra_features)):
                for i in range(normal_L):
                    f_pos[:, i] *= f_pos[:, normal_L + j]
                    f_neg[:, i] *= f_neg[:, normal_L + j]
            return (f_pos[pos_inverse, :normal_L], f_neg[neg_inverse, :normal_L], pos_slice)
        else:
            return (f_pos[pos_inverse, :], f_neg[neg_inverse, :], pos_slice)

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
            #     cur.append(pos)
            #     c_labels.append(1)
            #     cur.append(neg)
            #     c_labels.append(-1)
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
                       row['cor'], row['offset']] + list(additional_info)
            print("\t".join(map(str, out_row)), file=file_desc)
