#!/usr/bin/python
# -------------------------------------------------------------
# Utility to perform incremental PCA on a matrix in hdf5 format
# Jun. 16th, 2020
# Carles Boix
# -------------------------------------------------------------
import re
import os
import gc
import sys
import gzip
import fire
import time
import pickle
import pandas as pd
import numpy as np
import h5py

from scipy.io import mmread, mminfo
import scipy.sparse as sparse
from sklearn.utils import sparsefuncs
from sklearn.decomposition import IncrementalPCA

# For dim reduction:
from sklearn.manifold import TSNE
import umap

import socket

domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')
# from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from matplotlib import gridspec
import seaborn as sns


# Try incremental PCA with sklearn:
# lddir = 'linking_data/'
# datasuf = '_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged'
# fullpref = lddir + 'merged6mark' + datasuf
# merge_hdf5 = fullpref + '.hdf5'
# h5ad_file = fullpref + "_scmatrix.h5ad"
# prefix = fullpref + "_scrun_"
# k = 30
# merge_xpca = fullpref + '_pca_k' + str(k) + '.tsv.gz'

# h5 = h5py.File(merge_hdf5, 'r')
# data = h5['alldata']
# n = data.shape[0] # how many rows we have in the dataset
# chunk_size = 25000 # how many rows we feed to IPCA at a time, the divisor of n

class run_pca(object):
    def __init__(self, hdf5file, outprefix, k=10,
                 chunk_size=25000, batch_size=1000,
                 precomputed=False, plotpcs=False, transpose=False,
                 tsne=False, umap=False, normalize=False):
        self.hdf5file = str(hdf5file)
        self.outprefix = str(outprefix)
        self.k = int(k)
        self.chunk_size = int(chunk_size)
        self.batch_size = int(batch_size)
        # After running PCA:
        self.precomputed = precomputed
        self.plotpcs = str(plotpcs)
        self.normalize = normalize
        self.transpose = transpose
        # Run tsne, umap or not:
        self.tsne = tsne
        self.umap = umap
        # Files for PCA:
        self.pcafile = self.outprefix + '.pca_k' + str(self.k) + '.tsv.gz'

    # Main workflow:
    def main(self):
        if self.precomputed:
            self.read_results()
        else:
            self.perform_pca()
            self.write_results()
        if self.plotpcs:
            self.plot_pcs(filename=self.outprefix + '_pcs.png', k=self.k)
        if self.umap:
            self.run_umap()
            self.plot_dimred(self.umap_U, filename=self.outprefix + '_umap.png',
                             redtype='UMAP')
        if self.tsne:
            self.run_tsne()
            self.plot_dimred(self.tsne_U, filename=self.outprefix + '_tsne.png',
                             redtype='t-SNE')

    def perform_pca(self):
        self.h5 = h5py.File(self.hdf5file, 'r')
        self.data = self.h5['matrix']
        self.n = self.data.shape[0]
        # NOTE: Batch should be bigger (5x?) than the number of features
        self.niter = self.n // self.chunk_size
        self.ipca = IncrementalPCA(n_components=self.k,
                                   batch_size=self.batch_size)
        print("[STATUS] Running incremental fit:")
        t0 = time.time()
        for i in range(0, self.niter):
            t1 = time.time() - t0
            print(str(i) + ' (' +  str(round(i / self.niter *100, 1)) +
                  '%) at' + str(round(t1, 2)) +  's')
            i0 = i * self.chunk_size
            i1 = (i + 1) * self.chunk_size
            self.ipca.partial_fit(self.data[i0:i1])
        # Apply the PCA decomp:
        self.X_pca = np.zeros((self.n, self.k))
        print("[STATUS] Applying fit:")
        t0 = time.time()
        for i in range(0, self.niter):
            t1 = time.time() - t0
            print(str(i) + ' (' +  str(round(i / self.niter *100, 1)) +
                  '%) at' + str(round(t1, 2)) +  's')
            i0 = i * self.chunk_size
            i1 = (i + 1) * self.chunk_size
            self.X_pca[i0:i1,:] = self.ipca.fit_transform(self.data[i0:i1])
        print("[STATUS] Writing PCs:")

    def write_results(self):
        print("[STATUS] Writing out results")
        gzipped_write(self.X_pca, self.pcafile)

    def read_results(self):
        print("[STATUS] Reading in results")
        self.X_pca = pd.read_csv(self.pcafile, header=None, sep="\t").to_numpy()

    def run_umap(self):
        t1 = time.time()
        print("[STATUS] Performing UMAP")
        fit = umap.UMAP()
        self.umap_U = fit.fit_transform(self.X_pca)
        print('[STATUS] Finished in ' + str(round(time.time() - t1, 2)) + "s")
        print("[STATUS] Writing out UMAP results")
        umapfile = self.outprefix + '.umap.tsv.gz'
        gzipped_write(self.umap_U, umapfile)

    def run_tsne(self):
        t1 = time.time()
        print("[STATUS] Performing t-SNE")
        fit = TSNE(n_components=2)
        self.tsne_U = fit.fit_transform(self.X_pca)
        print('[STATUS] Finished in ' + str(round(time.time() - t1, 2)) + "s")
        print("[STATUS] Writing out t-SNE results")
        tsnefile = self.outprefix + '.tsne.tsv.gz'
        gzipped_write(self.tsne_U, tsnefile)

    def plot_dimred(self, mat, filename, redtype=None, title=None):
        sns.set(font_scale=1.1)
        # Plot current trace:
        fig = plt.figure(figsize=(8, 8))
        ax = plt.gca()
        # Plot % change:
        if title is not None:
            ax.set_title(title)
        ax.set_facecolor('white')
        plt.scatter(mat[:,0], mat[:,1], s=.1)
        if redtype is not None:
            plt.ylabel(redtype + ' 1')
            plt.xlabel(redtype + ' 2')
        plt.tight_layout()
        fig = plt.gcf()
        fig.savefig(filename, dpi=350, bbox_inches='tight')

    # Plot main (all) principal components
    def plot_pcs(self, filename, k):
        c = 5
        r = k // c
        sns.set(font_scale=1.1)
        fig = plt.figure(figsize=(2 * c, 2 * r))
        gs = gridspec.GridSpec(r, c)
        for p in range(5):
            ax = plt.subplot(gs[p])
            ax.set_facecolor('white')
            if ((2 * p + 1) > self.k):
                self.plot_pcpair(2 * p, 2 * p + 1, ax)
        plt.tight_layout()
        fig = plt.gcf()
        fig.savefig(filename, dpi=350, bbox_inches='tight')

    def plot_pcpair(self, i, j, ax):
        ax.set_facecolor('white')
        plt.scatter(self.X_pca[:,i], self.X_pca[:,j], s=.1)
        plt.ylabel('PC ' + str(i + 1))
        plt.xlabel('PC ' + str(j + 2))


def gzipped_write(matrix, filename):
    with gzip.open(filename, 'wb') as f:
        preformatted_write(matrix, f, encode=True)

# Faster than for loop or np.savetxt:
def preformatted_write(mat, f, fmtstring=None, encode=False):
    if fmtstring is None:
        if len(mat.shape) == 1:
            fmtstring = '%g'
        else:
            fmtstring = '\t'.join(['%g']*mat.shape[1])
    fmt = '\n'.join([fmtstring]*mat.shape[0])
    data = fmt % tuple(mat.ravel())
    if encode:
        data = data.encode('utf-8')
    f.write(data)


if __name__ == "__main__":
    fire.Fire(run_pca)
