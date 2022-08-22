#!/usr/bin/python
# ----------------------------------
# Auxiliary functions for processing
# ----------------------------------
import os
import gzip
import numpy as np
import six.moves.cPickle as pickle
from scipy.sparse import csr_matrix

# Clustering:
from scipy.cluster import hierarchy as sch
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

rng = np.random


def diagonalize_mat(mat, cutoff=0.05):
    Z = mat > cutoff
    colsums = np.sum(Z, axis=1)
    cols = np.array([i for i, v in enumerate(colsums) if v > 0])
    maxids = [np.max(np.where(Z[i, :])[0]) for i in cols]
    diag_idx = np.argsort(maxids)
    reord = mat[cols, :]
    reord = reord[diag_idx, ]
    return(reord, diag_idx, cols)


# TODO fixme better
def pairwise_jaccard(X):
    # X = X.astype(bool).astype(int)
    intrsct = X.dot(X.T)
    # row_sums = np.squeeze(np.array(np.sum(X, axis=1)))
    row_sums = intrsct.diagonal()
    unions = row_sums[:, None] + row_sums - intrsct
    dist = 1.0 - intrsct / unions
    return(dist)


# Metric: euclidean or cosine works ok. Ward is best linkage method.
def reorder_by_linkage(mat, method='ward', metric='cosine', msg=False):
    if msg:
        print("Reordering inner dim by linkage...")
    mat = mat + np.abs(0.0001 * rng.normal(0, 1, size=mat.shape))
    D = pdist(mat, metric=metric)
    Z = linkage(D, method=method)
    # Z = optimal_leaf_ordering(Z, D)
    row_idx = sch.leaves_list(Z)
    if msg:
        print("Reordering columns by linkage...")
    D = pdist(mat.T, metric=metric)
    Z = linkage(D, method=method)
    # Z = optimal_leaf_ordering(Z, D)
    col_idx = sch.leaves_list(Z)
    reord = mat[row_idx[:, np.newaxis], col_idx]
    return(reord, row_idx, col_idx)


# TODO: Write bed file as gzipped?
# TODO: Write merged and separate (awk?)
def write_bedfile(filename, chrom, ind, clust, width=200):
    # TODO: REMOVE THE PREVIOUS WIDTH!
    if len(ind) == len(clust):
        with open(filename, 'wb') as outfile:
            for i in range(len(ind)):
                line = [chrom[i],
                        str(ind[i] * width + 1),
                        str((ind[i] + 1) * width),
                        str(clust[i])]
                line = "\t".join(line) + "\n"
                outfile.write(line)
                pass
    else:
        print("Index and clusters not same length")


def load_pickle(filename):
    with open(filename, 'rb') as infile:
        matrix = pickle.load(infile)
    return matrix


def load_pickle_gzip(filename):
    with gzip.open(filename, 'rb') as infile:
        matrix = pickle.load(infile)  # , encoding='latin1')
    return matrix


def load_pickle_gzip_latin(filename):
    with gzip.open(filename, 'rb') as infile:
        matrix = pickle.load(infile, encoding='latin1')
    return matrix


# TODO add formatter sprintf
def write_list(filename, ll):
    with open(filename, 'w') as outfile:
        for item in ll:
            outfile.write(str(item) + "\n")


def save_pickle(filename, matrix):
    try:
        with open(filename, 'wb') as outfile:
            pickle.dump(matrix, outfile, protocol=2)
    except:
        with open(filename, 'wb') as outfile:
            pickle.dump(matrix, outfile, protocol=4)


def save_pickle_gzip(filename, matrix):
    try:
        with gzip.open(filename, 'wb') as outfile:
            pickle.dump(matrix, outfile, protocol=2)
    except:
        with gzip.open(filename, 'wb') as outfile:
            pickle.dump(matrix, outfile, protocol=4)


def save_sparse_csr(filename, array):
    np.savez(filename, data=array.data, indices=array.indices,
             indptr=array.indptr, shape=array.shape)


def load_sparse_csr(filename):
    loader = np.load(filename + '.npz')
    return csr_matrix((loader['data'], loader['indices'], loader['indptr']),
                      shape=loader['shape'])


def load_file_save_sparse(path):
    sparse_file = path + "_csr.cp"
    if os.path.isfile(sparse_file):
        with open(sparse_file, "rb") as fp:
            X = pickle.load(fp)
    elif os.path.isfile(sparse_file + ".gz"):
        with gzip.open(sparse_file + ".gz", "rb") as fp:
            X = pickle.load(fp)
    else:
        X = np.load(file=path + ".npy")
        X = csr_matrix(X)
        save_pickle(sparse_file, X)
    return(X)


# Write a list into a file:
def savelist(x, fname):
    with open(fname, 'w') as f:
        for item in x:
            f.write(str(item) + '\n')


# Return the averaged windows (200bp):
def read_bw_avg_idlist(file, idlist, window=8):
    # NOTE: idlist must be sorted/unique:
    if not is_sorted(idlist):
        idlist = np.unique(np.sort(idlist))
    if (len(idlist) % window) != 0:
        print("[ERROR] IDLIST MUST BE MULTIPLE OF WINDOWSIZE")
    array = []
    with gzip.open(file, 'rt') as f:
        # Read/skip two header lines:
        header1 = next(f)
        header2 = next(f)
        ind = 0
        widx = 0
        val = 0
        bufsize = 65536
        while True:
            lines = f.readlines(bufsize)
            if not lines or len(idlist) == 0:
                break
            for line in lines:
                # Add only if in list, othw. keep going:
                if ind == idlist[0]:
                    line = line.strip("\n")
                    val = val + float(line)
                    widx = widx + 1  # Number of values added
                    # Average every eight values:
                    if widx == window:
                        array.append(val / window)
                        val = 0
                        widx = 0
                    idlist = idlist[1:]
                ind = ind + 1
                # Stop if no more ids to pull
                if len(idlist) == 0:
                    break
    narray = np.array(array)
    return([narray, header1, header2])


# Lazy/fast check of sortedness
def is_sorted(a):
    for i in range(a.size-1):
        if a[i+1] < a[i]:
            return(False)
    return(True)
