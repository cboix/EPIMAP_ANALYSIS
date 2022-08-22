#!/usr/bin/python
# --------------------------------------------------------------------
# Script to return genome-wide statistics on variabilty of hdf5 matrix
# --------------------------------------------------------------------
import fire
import h5py
import gzip
import numpy as np
import pandas as pd


class calculate_mat_stats(object):
    def __init__(self, hdf5file, outprefix):
        self.hdf5file = hdf5file
        self.outprefix = outprefix

    def main(self):
        self.read_matrix()
        self.calculate_statistics()
        self.write_statistics()

    def read_matrix(self):
        with h5py.File(self.hdf5file, 'r') as f:
            self.X = f['matrix'][:,:]
        print('[STATUS] Read in matrix of shape:', self.X.shape)

    def calculate_statistics(self):
        self.dd = {}
        self.dd['Mean'] = np.mean(self.X, axis=1)  # Mean
        self.dd['Max'] = np.max(self.X, axis=1)  # Mean
        self.dd['Std'] = np.std(self.X, axis=1)  # Std. dev
        self.dd['CoV'] = self.dd['Std'] / self.dd['Mean']  # Coefficient of variation
        self.dd['nsamp'] = np.sum(self.X >= 2, axis=1)  # Num samples per peak
        # Turn into dataframe:
        self.df = pd.DataFrame(self.dd)

    def write_statistics(self):
        self.df.to_csv(self.outprefix + "_df.tsv.gz",
                       index=False, sep="\t")

if __name__ == "__main__":
    fire.Fire(calculate_mat_stats)

