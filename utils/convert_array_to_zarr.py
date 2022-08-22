#!/usr/bin/python
# -----------------------------------------------------------
# Utility to save an array as a zarr file
# To support: mtx, cp, and txt/tsv formats
# April 12, 2020
# Carles Boix
# -----------------------------------------------------------
import os
import re
import gzip
import fire # CLI
import time
# Compression schemes:
from numcodecs import Blosc

# Data formats:
import zarr
import h5py
import pickle
import numpy as np
import pandas as pd
import scipy.sparse as sparse
from scipy.io import mminfo


class convert_to_zarr(object):
    def __init__(self, infile, outfile, chunk1=None, chunk2=None,
                 auto_chunk=False, clevel=4,
                 make_zarr=True, make_h5py=False):
        self.infile = str(infile)
        self.outfile = str(outfile)
        self.clevel = clevel
        self.chunk1 = chunk1
        self.chunk2 = chunk2
        self.auto_chunk = auto_chunk
        self.make_zarr = make_zarr
        self.make_h5py = make_h5py

    # Main workflow:
    def main(self):
        self.load_file()
        if not self.auto_chunk:
            self.determine_chunksize()
        if self.make_zarr:
            self.save_zarr()
        if self.make_h5py:
            self.save_h5py()


    def load_file(self):
        t1 = time.time()
        if re.search('\\.mtx', self.infile):
            self.read_mtx()
        elif re.search('\\.cp', self.infile):
            self.read_cp()
        elif re.search('\\.tsv', self.infile) or re.search(
            '\\.txt', self.infile):
            self.read_tsv()
        else:
            print("[WARNING] FILETYPE NOT RECOGNIZED DEFAULTING TO CSV")
            self.filetype = 'txt'
        print("[STATUS] Read in mtx file with dimensions:", self.X.shape)
        print('[STATUS] Finished in ' + str(round(time.time() - t1, 2))+"s")

    # Read in file (or pickle)
    def read_mtx(self):
        t1 = time.time()
        print("[STATUS] Reading in mtx file:", self.infile)
        self.mtxinfo = mminfo(self.infile)
        print("[STATUS] Matrix info:", self.mtxinfo)
        self.xdf = pd.read_csv(self.infile, skiprows=3, header=None, sep=" ",
                               names=['i','j','v']).to_numpy()
        print(self.xdf[0:5,:])
        sparsity = self.mtxinfo[2] / (self.mtxinfo[0] * self.mtxinfo[1])
        print('sparsity', sparsity)
        self.X = sparse.coo_matrix((self.xdf[:,2],
                                    (self.xdf[:,0] - 1, self.xdf[:,1] - 1)),
                                   shape=(self.mtxinfo[0], self.mtxinfo[1]))
        # Only keep as sparse if at least 50% sparse:
        if sparsity > 0.50:
            self.X = self.X.toarray()

    def read_tsv(self):
        print("[STATUS] Loading from tsv with pd.read_csv")
        self.X = pd.read_csv(self.infile, header=None, sep="\t").to_numpy()

    def read_cp(self):
        print("[STATUS] Loading from pickle")
        self.X = load_pickle_gzip_latin(self.infile)

    def determine_chunksize(self):
        print("Choosing chunksize")
        if self.chunk1 is None:
            if self.X.shape[0] < 1000:
                self.chunk1 = self.X.shape[0]
            elif self.X.shape[0] < 10000:
                self.chunk2 = 1000
            else:
                self.chunk2 = 10000
        if self.chunk2 is None:
            if self.X.shape[1] < 1000:
                self.chunk2 = self.X.shape[1]
            elif self.X.shape[1] < 10000:
                self.chunk2 = 1000
            else:
                self.chunk2 = 10000
        print(self.chunk1, self.chunk2)

    def save_zarr(self):
        print("[STATUS] Saving as zarr:", self.outfile)
        compressor = Blosc(cname='zstd', clevel=self.clevel, shuffle=Blosc.BITSHUFFLE)
        if type(self.X) == sparse.csr_matrix:
            # TODO: FIXME
            print("SPARSE NOT IMPLEMENTED YET")
            a = sparse.COO(coords, data, shape=shape)
            z = zarr.zeros(shape=a.shape, chunks=a.ndim*(1,),
                           dtype=a.dtype, compressor=compressor)
            z.set_coordinate_selection(tuple(a.coords), a.data)
        else:
            if self.chunk1 is None and self.chunk2 is None:
                self.Z = zarr.array(self.X, store=self.outfile,
                                    overwrite=True, compressor=compressor)
            else:
                self.Z = zarr.array(self.X, chunks=(self.chunk1, self.chunk2),
                                    store=self.outfile, overwrite=True,
                                    compressor=compressor)
        print(self.Z.info)

    def save_h5py(self):
        print("Converting to dense", self.outfile)
        if type(self.X) == sparse.csr_matrix:
            print("From CSR")
            self.X = self.X.toarray()
            print('sum', np.sum(self.X))
            print('nnz', np.sum(self.X > 0))
        elif type(self.X) == sparse.coo_matrix:
            print("From COO")
            self.X = self.X.toarray()
            print('sum', np.sum(self.X))
            print('nnz', np.sum(self.X > 0))
        print("[STATUS] Saving as hdf5:", self.outfile)
        self.hf = h5py.File(self.outfile, 'w')
        self.hf.create_dataset('matrix', data=self.X,
                compression='gzip', compression_opts=self.clevel)
        self.hf.close()


def load_pickle_gzip_latin(filename):
    with gzip.open(filename, 'rb') as infile:
        matrix = pickle.load(infile, encoding='latin1')
    return(matrix)

if __name__ == "__main__":
    fire.Fire(convert_to_zarr)
