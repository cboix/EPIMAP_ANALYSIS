#!/usr/bin/env python3
# Author: Benjamin T. James
import sys
import h5py
import numpy as np

class mark:
    def __init__(self, mark_list, sample_list_indices, mark_neg_idx_list=None, log_file=sys.stderr):
        self.data = None
        self.log = log_file
        if type(sample_list_indices) != list or (len(sample_list_indices) > 0 and type(sample_list_indices[0]) is not int):
            raise ValueError('%s is not a valid list of integers', sample_list_indices)
        if mark_list is None:
            mark_list = []
        for i, f in enumerate(mark_list):
            print("Reading in", f, file=self.log)
            fd = h5py.File(f, 'r')
            if self.data is None:
                self.data = np.zeros((len(mark_list), fd['matrix'].shape[0]))
                self.ndata = np.zeros((len(mark_list), fd['matrix'].shape[0]))
            self.data[i, :] = np.mean(fd['matrix'][:, sample_list_indices], axis=1)


            if mark_neg_idx_list is not None:
                ### if provided, use negative sample instead of shuffling indices
                self.ndata[i, :] = np.mean(fd['matrix'][:, mark_neg_idx_list], axis=1)
            else:
                self.ndata[i, :] = self.data[i, :]
                np.random.shuffle(self.ndata[i, :])
            #self.ndata[i, :] = np.median(fd['matrix'][:, ])
            fd.close()
        if self.data is None:
            self.data = np.zeros((0, 0))

    def featurize(self, key='mind'):
        out = []
        for i in range(self.data.shape[0]):
            fn = lambda p, n: (self.data[i, p[key].values],
                               self.ndata[i, p[key].values])
            # fn = lambda p, n: (self.data[i, p[key].values],
            #                    self.data[i, np.random.randint(0, self.data.shape[1], len(p[key].values))])
            out.append(fn)
        return out
