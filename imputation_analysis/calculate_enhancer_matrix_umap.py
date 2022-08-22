# !/usr/bin/python
# ----------------------------------------------
# Use scanpy to perform dimensionality reduction
# and  clustering on enhancer matrix
# ----------------------------------------------
import glob
import h5py
from scipy import sparse
import re
import os
from sklearn.utils import sparsefuncs
import gzip
import pickle
import numpy as np
import pandas as pd
import time
import gc
import fbpca
import sys
import h5py
from sklearn.decomposition import IncrementalPCA

# For KNN + conn
import umap
from umap.umap_ import nearest_neighbors
from sklearn.utils import check_random_state
from types import MappingProxyType
from umap.umap_ import fuzzy_simplicial_set

# For clustering:
import igraph as ig
import leidenalg
import hdbscan

# For plotting:
import socket
domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')

from matplotlib import pyplot as plt
import seaborn as sns

import scanpy.api as sc
import anndata
from scanpy.preprocessing import normalize_total

# Scanpy plotting settings:
sc.settings.verbosity = 2  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=300)  # set sufficiently high resolution for saving



# ----------
# Functions:
# ----------
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

def gzipped_write(matrix, filename):
    with gzip.open(filename, 'wb') as f:
        preformatted_write(matrix, f, encode=True)

def plot_dimred(mat, filename, redtype=None, title=None,
                values=None, cmap=None, size=.1):
    sns.set(font_scale=1.1)
    # Plot current trace:
    fig = plt.figure(figsize=(8, 8))
    ax = plt.gca()
    # Plot % change:
    if title is not None:
        ax.set_title(title)
    ax.set_facecolor('white')
    if values is None:
        plt.scatter(mat[:,0], mat[:,1], s=size)
    else:
        plt.scatter(mat[:,0], mat[:,1], s=size,
                    c=values, cmap=cmap)
    if redtype is not None:
        plt.ylabel(redtype + ' 1')
        plt.xlabel(redtype + ' 2')
    plt.tight_layout()
    fig = plt.gcf()
    fig.savefig(filename, dpi=350, bbox_inches='tight')

def load_pickle_gzip_latin(filename):
    with gzip.open(filename, 'rb') as infile:
        matrix = pickle.load(infile, encoding='latin1')
    return(matrix)

# For labels write:
def aggregate_labels_over_matrix(matrix, labels, hdb=False):
    NLAB = np.max(labels) + 1 + 1 * hdb
    avgmat = np.zeros((NLAB, matrix.shape[1]))
    for i in range(NLAB):
        print(i)
        ind = np.where(labels == (i - 1 * hdb))[0]
        avgmat[i,:] = np.array(np.mean(matrix[ind,:], axis=0)[0])[0]
    return(avgmat)

def calc_write_avgmatrix(lblset, prefix, fullmat, hdb=False):
    groupfile = prefix + "." + lblset + '.tsv.gz'
    avgfile = prefix + "." + lblset + '.avg.tsv.gz'
    labels = pd.read_csv(groupfile, header=None, sep="\t").to_numpy().T[0]
    avgmat = aggregate_labels_over_matrix(fullmat, labels, hdb=hdb)
    gzipped_write(avgmat, avgfile)



# ------------
# Data import:
# ------------
# OPTIONS:
matdir = 'masterlist_matindices/'
matpref = matdir + "matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH"
indfile = matpref + "_masterlist_indices.cp.gz"
use_enh = False
if use_enh:
    matfile = matpref + "_matrix_csr.cp.gz"
    attrfile = matpref + "_matrix_attr.cp.gz"
    h5ad_file = matpref + "_matrix.h5ad"
    prefix = matpref + "_scrun"
else:
    acpref = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged'
    ac_hd5 = acpref + '.hdf5'
    acnames = 'linking_data/mark_matrix_names.txt'
    h5ad_file = acpref + "_scmatrix.h5ad"
    prefix = acpref + "_scrun"


# Instead read in the full matrix:
if not os.path.isfile(h5ad_file):
    enhind = load_pickle_gzip_latin(indfile)
    if use_enh:
        # Load in the data:
        attr = load_pickle_gzip_latin(attrfile)
        X = load_pickle_gzip_latin(matfile)
        names = attr['names']
    else:
        names = pd.read_csv(acnames, header=None).to_numpy().T[0]
        with h5py.File(ac_hd5, 'r') as f:
            X = f['matrix'][:,:]
    print(X.shape)
    X = X[enhind,:]
    print(X.shape)
    # Make the adata object:
    d = {}
    d['obs'] = ['d' + str(i) for i in range(X.shape[0])]
    d['var'] = names
    d['X'] = X
    adata = anndata.AnnData(**d)
    adata.obs_names = d['obs']
    adata.var_names = d['var']
    adata.var.columns = ['varnames']
    adata.obs.columns = ['obsnames']
    del(X)
    gc.collect()
    # Write the adata object:
    adata.write(h5ad_file)
else:
    adata = sc.read_h5ad(h5ad_file)


# ----------------------------
# On multi-mark data, no norm:
# ----------------------------
# Try incremental PCA with sklearn:
lddir = 'linking_data/'
datasuf = '_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged'
fullpref = lddir + 'merged6mark' + datasuf
merge_hdf5 = fullpref + '.hdf5'
h5ad_file = fullpref + "_scmatrix.h5ad"
prefix = fullpref + "_scrun_"
k = 30
merge_xpca = fullpref + '_pca_k' + str(k) + '.tsv.gz'

if not os.path.isfile(h5ad_file):
    if not os.path.isfile(merge_xpca):
        h5 = h5py.File(merge_hdf5, 'r')
        data = h5['alldata']
        n = data.shape[0] # how many rows we have in the dataset
        chunk_size = 25000 # how many rows we feed to IPCA at a time, the divisor of n
        # NOTE: Batch should be bigger (5x?) than the number of features
        # Storage takes a while, 180 Takes about 6 hrs (2-??pm ish - with 1000=batch and 3.6M * 5K and K=30)
        ipca = IncrementalPCA(n_components=k, batch_size=1000)
        for i in range(0, n//chunk_size):
            print(i)
            ipca.partial_fit(data[i*chunk_size : (i+1)*chunk_size])
        # Apply the PCA decomp:
        X_pca = np.zeros((n, k))
        for i in range(0, n//chunk_size):
            print(i)
            X_pca[(i*chunk_size):((i+1)*chunk_size),:] = ipca.fit_transform(
                data[(i*chunk_size):((i+1)*chunk_size)])
        # SAVE:
        gzipped_write(X_pca, merge_xpca)
    else:
        X_pca = pd.read_csv(merge_xpca, header=None, sep="\t").to_numpy()
    d = {}
    d['obs'] = ['d' + str(i) for i in range(X_pca.shape[0])]
    d['var'] = ['pc' + str(i) for i in range(X_pca.shape[1])]
    d['X'] = X_pca
    adata = anndata.AnnData(**d)
    adata.obs_names = d['obs']
    adata.var_names = d['var']
    adata.var.columns = ['varnames']
    adata.obs.columns = ['obsnames']
    adata.obsm['X_pca'] = X_pca
    gc.collect()
    # Write the adata object:
    adata.write(h5ad_file)
else:
    adata = sc.read_h5ad(h5ad_file)


# --------------
# Preprocessing:
# --------------
if not hasattr(adata, 'obsm'):
    # normalize_total(adata) # Destroys signal
    # sc.pp.log1p(adata) # Unclear if necessary
    # PCA:
    # Timing: 3 minutes with no zero_center? norm is affecting.
    t1 = time.time()
    k = 50
    pcafile = prefix + '.X_pca.tsv.gz'
    sc.pp.pca(adata, n_comps=k, zero_center=False)
    gzipped_write(adata.obsm['X_pca'], pcafile)
    print('PCA',time.time() - t1)

# NOTE: NN=50 on enh ~ 45min
# NN=200 = 2:20 on ac
t1 = time.time()
# sc.pp.neighbors(adata) # Default of 15 may be too few
NN = 200 # Dense 3M * 200 doesn't seem to fit in memory
NN = 20
sc.pp.neighbors(adata, n_neighbors=NN)
print('Neighbors',time.time() - t1)
adata.write(h5ad_file)

niter = 50
RES = 1
prefstr = '_r' + str(RES) + '_n' + str(niter)
sc.tl.leiden(adata, n_iterations=niter, resolution=RES)
adata.write(h5ad_file)
labels = adata.obs['leiden'].to_numpy().astype(int)
gzipped_write(labels, prefix + '.leiden' + prefstr + '.tsv.gz')
print('Leiden',time.time() - t1)

xp = adata.obsm['X_pca']

sns.set(font_scale=1.1)
# Plot current trace:
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
ax.set_facecolor('white')
plt.scatter(xp[:,0], xp[:,1], s=.1)
plt.ylabel('PC 1')
plt.xlabel('PC 2')
plt.tight_layout()
fig = plt.gcf()
fig.savefig(prefix + "pc_12.png", dpi=350, bbox_inches='tight')

sns.set(font_scale=1.1)
# Plot current trace:
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
ax.set_facecolor('white')
plt.scatter(xp[:,2], xp[:,3], s=.1)
plt.ylabel('PC 1')
plt.xlabel('PC 2')
plt.tight_layout()
fig = plt.gcf()
fig.savefig(prefix + "pc_34.png", dpi=350, bbox_inches='tight')

sns.set(font_scale=1.1)
# Plot current trace:
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
ax.set_facecolor('white')
plt.scatter(xp[:,3], xp[:,4], s=.1)
plt.ylabel('PC 1')
plt.xlabel('PC 2')
plt.tight_layout()
fig = plt.gcf()
fig.savefig(prefix + "pc_56.png", dpi=350, bbox_inches='tight')


sns.set(font_scale=1.1)
# Plot current trace:
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
ax.set_facecolor('white')
plt.scatter(xp[:,8], xp[:,9], s=.1)
plt.ylabel('PC 1')
plt.xlabel('PC 2')
plt.tight_layout()
fig = plt.gcf()
fig.savefig(prefix + "pc_90.png", dpi=350, bbox_inches='tight')



sns.set(font_scale=1.1)
# Plot current trace:
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
ax.set_facecolor('white')
plt.scatter(xp[:,18], xp[:,19], s=.1)
plt.ylabel('PC 1')
plt.xlabel('PC 2')
plt.tight_layout()
fig = plt.gcf()
fig.savefig(prefix + "pc_1920.png", dpi=350, bbox_inches='tight')

sns.set(font_scale=1.1)
# Plot current trace:
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
ax.set_facecolor('white')
plt.scatter(xp[:,48], xp[:,49], s=.1)
plt.ylabel('PC 1')
plt.xlabel('PC 2')
plt.tight_layout()
fig = plt.gcf()
fig.savefig(prefix + "pc_4950.png", dpi=350, bbox_inches='tight')



# Spectral doesn't work --> too many separate connected components?
# sc.tl.umap(adata, maxiter=1, init_pos='spectral') # takes forever and crashes.
# sc.tl.umap(adata, maxiter=1, init_pos='random') # Takes 4:37 min, random square.
# sc.tl.umap(adata, maxiter=2, init_pos='random') # Takes 4:38 min??
# sc.tl.umap(adata, maxiter=10, init_pos='random') # Takes 6:27
# sc.tl.umap(adata, init_pos='random') # No iter limit. Takes XXX
sc.tl.umap(adata, init_pos='spectral') # No iter limit. Takes XXX


# UMAP parameters:
# adata, min_dist=0.5, spread=1.0, n_components=2, maxiter=None,
# alpha=1.0, gamma=1.0, negative_sample_rate=5, init_pos='spectral',
# random_state=0, a=None, b=None, copy=False, method='umap', neighbors_key=None
imgpref = re.sub(".*/","", prefix)
sc.pl.umap(adata, frameon=False, save=imgpref + '_umap.png')

sc.pl.umap(adata, color='leiden', legend_loc='on data',
           frameon=False, save=imgpref + '_leiden.png')

adata.write(h5ad_file)
gzipped_write(adata.obsm['X_umap'], prefix + '.umap.tsv.gz')
print('UMAP',time.time() - t1)

# sc.tl.rank_genes_groups(adata, 'leiden')
# sc.pl.rank_genes_groups(adata, save='.pdf') # Plot
# adata.write(h5ad_file)
# print('Ranking',time.time() - t1)

# Get the average values for these clusters:
lblset = 'leiden' + prefstr
calc_write_avgmatrix(lblset, prefix=prefix, fullmat=adata.X)







