#!/usr/bin/python
# ---------------------------------------------
# Standalone utility to turn df into sparse csr
# Format of dataframe:
# col1 col2 value
# Tool returns attributes and sparse matrix
# ---------------------------------------------
import os
import fire
import csv
import gzip
import numpy as np
from scipy.sparse import csr_matrix, hstack, vstack
# import argparse
from tqdm import tqdm

# Load config into os.environ:
# import parse_config
import auxfunc_masterlist as AUX


class txtdf_to_sparse(object):
    def __init__(self, txtfile, outprefix, binary=False, verbose=False):
        # Arguments:
        self.txtfile = txtfile
        self.outprefix = outprefix
        self.binary = binary
        # https://stackoverflow.com/questions/5980042/
        if self.verbose:
            def verboseprint(self, *args):
                # Print each argument separately so caller doesn't need to
                # stuff everything to be printed into a single string
                for arg in args:
                    print(arg)
                print
        else:
            self.verboseprint = lambda *a: None  # do-nothing function

    def main(self):
        self.load_file()
        self.convert_indices()
        self.make_matrix()
        self.write_matrix_attr()

    def load_file(self):
        self.icol = []
        self.jcol = []
        self.vals = []
        with open(self.txtfile, 'r') as f:
            for line in f:
                row = line.rstrip("\n").split("\t")
                self.icol.append(row[0])
                self.jcol.append(row[1])
                # TODO: Specify dtype
                self.vals.append(float(row[2]))
        print("[STATUS] Read in txtfile: " + str(len(self.vals)) + " values.")

    def convert_indices(self):
        self.ui = np.unique(self.icol)
        self.uj = np.unique(self.jcol)

    def make_matrix(self):
        pass

    def write_matrix_attr(self):
        pass

if __name__ == "__main__":
    fire.Fire(txtdf_to_sparse)
