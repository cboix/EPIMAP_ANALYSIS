#!/usr/bin/python
# ----------------------------------------
# Script to compress all BW files
# for a specific mark # an h5py dataframe
# ----------------------------------------
import os
import socket
import fire
import glob
import csv
import gzip
import re
import numpy as np
import h5py
import time
import six.moves.cPickle as pickle
from scipy.sparse import csr_matrix  # , hstack, vstack
# import argparse
from tqdm import tqdm


class compress_bw(object):
    def __init__(self, mark, outdir='hdf5_bw/', dataset='imputed', chrom=None):
        # Directories:
        self.domain = socket.getfqdn()
        if 'broadinstitute.org' in self.domain:
            self.HOME = '/broad/compbio/cboix'
        else:
            self.HOME = os.getenv("HOME")
        self.datadir = self.HOME + '/EPIMAP_ANALYSIS/db/ChromImpute/'
        # Arguments:
        self.mark = mark
        self.outdir = self.datadir + outdir
        self.dataset = dataset
        # Either process one chrom or all:
        self.chrom = chrom
        if self.chrom is None:
            self.chrlist = ['chr' + str(i+1) for i in range(22)] + ['chrX']
        else:
            if type(chrom) == str:
                self.chrlist = [chrom]
            else:
                self.chrlist = chrom
                # Directory of files:
        if self.dataset == 'imputed':
            self.bwdir = self.datadir + 'imputed/'
            self.pattern = self.chrlist[0] + "_impute" + "_*" + self.mark + "*.wig.gz"
        elif self.dataset == 'observed':
            self.bwdir = self.datadir + 'converted/'
            self.pattern = self.chrlist[0] + "_FINAL" + "_*.wig.gz"
        else:
            raise Exception("[ERROR] Do not recognize dataset: " +
                            str(self.dataset))
        self.prefix = self.dataset + '_' + str(self.mark)
        self.start = time.time()

    def main(self):
        self.setup()
        self.process()

    def setup(self):
        print("[STATUS] Searching directory for files")
        # Get list of files in directory:
        self.filelist = glob.glob(self.bwdir + self.pattern)
        self.basenames = [re.sub(".*/", "", filename)
                          for filename in self.filelist]
        self.suffixes = [re.sub(self.chrlist[0] + "_", "", basename)
                         for basename in self.basenames]
        self.ids = [re.sub(".wig.gz", "", suffix) for suffix in self.suffixes]
        self.ids = np.sort(np.unique(self.ids))
        print("[STATUS] There are " + str(len(self.ids)) +
              " prefixes of files matching " + self.pattern)

    def process(self):
        print("[STATUS] Processing files")
        # For each chromosome, join all files from the noted prefixes
        for chrom in self.chrlist:
            self.process_chrom(chrom)

    def process_concat_chrom(self, chrom):
        print("[STATUS] Processing " + chrom)
        j = 0
        for uid in tqdm(self.ids):
            filepref = self.bwdir + chrom + "_" + uid
            X = get_mat(filepref)
            if j == 0:
                self.FULL = X
                j = 1
            else:
                self.FULL = np.concatenate((self.FULL, X), axis=1)
        # Save array (in this case with h5py):
        self.save_chrom(chrom, self.FULL)

    def process_chrom(self, chrom):
        print("[STATUS] Processing " + chrom)
        j = 0
        for uid in tqdm(self.ids):
            filepref = self.bwdir + chrom + "_" + uid
            X = get_mat(filepref)  # Automatically creates CP file.

    # Export to hdf5 format:
    def save_chrom(self, chrom, X):
        identifier = chrom + "_" + self.mark + "_" + self.dataset
        print("[STATUS] Exporting " + identifier + " with shape " +
              str(X.shape) + " to HDF5 dataset.")
        self.f = h5py.File(self.outdir + identifier + '.hdf5', 'w')
        self.f.create_dataset(identifier, data=X, compression='gzip')
        self.f.close()


def get_mat(f):
    cp_file = f + '_csr.cp.gz'
    try:
        if os.path.isfile(cp_file):
            with gzip.open(cp_file, 'rb') as handle:
                X = pickle.load(handle)
        else:
            raise Exception("No such file: " + cp_file)
    except:
        with gzip.open(f + '.wig.gz', 'rt') as handle:
            reader = csv.reader(handle, delimiter='\t')
            line = next(reader)
            line = next(reader)
            X = []
            for line in reader:
                row = [float(x) for x in line]
                X.append(row)
            X = np.array(X)
        save_compress_pickle(cp_file, X)
    return(X)


def save_compress_pickle(filename, matrix):
    try:
        with gzip.open(filename, 'wb') as outfile:
            pickle.dump(matrix, outfile, protocol=2)
    except:
        with gzip.open(filename, 'wb') as outfile:
            pickle.dump(matrix, outfile, protocol=4)


def save_pickle(filename, matrix):
    try:
        with open(filename, 'wb') as outfile:
            pickle.dump(matrix, outfile, protocol=2)
    except:
        with open(filename, 'wb') as outfile:
            pickle.dump(matrix, outfile, protocol=4)


def save_sparse_csr(filename, array):
    np.savez(filename, data=array.data, indices=array.indices,
             indptr=array.indptr, shape=array.shape)


def load_sparse_csr(filename):
    loader = np.load(filename + '.npz')
    return csr_matrix((loader['data'], loader['indices'], loader['indptr']),
                      shape=loader['shape'])


def load_pickle(filename):
    with open(filename, 'rb') as infile:
        matrix = pickle.load(infile)
    return matrix


def load_file_save_sparse(path):
    sparse_file = path + "_csr.cp"
    if os.path.isfile(sparse_file):
        with open(sparse_file, "rb") as fp:
            X = pickle.load(fp)
    else:
        X = np.load(file=path + ".npy")
        X = csr_matrix(X)
        save_pickle(sparse_file, X)
    return(X)


if __name__ == "__main__":
    # Fire exposes vars and commands so we can run as:
    # python compress_bw_to_hd5py main --mark $MARK --dataset $DSET
    fire.Fire(compress_bw)
