#!/usr/bin/env python3
# Author: Benjamin T. James
# Adapted from code by Carles Boix
import sys
import re
import gzip
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import coo_matrix
import pickle
import os


# Load histone matrices to conform with others:
# prefix: Name up to extension
# name_order: intended order of the epigenomes
# binay: binarize matrix
def load_histonemat(prefix, name_order, binary=False):
    # Load names and get reordering:
    hnam = get_names(prefix)
    hnam = [re.sub(".*BSS", "BSS", n).split("_")[0] for n in hnam]
    name_order = hnam
    hord = np.array([np.where(np.array(hnam) == n)[0][0] for n in name_order])
    # Load and reorder matrix:
    tmpX = get_data(prefix)
    if binary:
        hXord = 1 * (tmpX[:, hord] > 0)
    else:
        hXord = tmpX[:, hord]
    del(tmpX)
    # Keep only the enhancers:
    # hXord = hXord[enhind,:]
    return(hXord)


def get_data(prefix):
    with gzip.open(prefix + "_csr.cp.gz", 'rb') as f:
        X = pickle.load(f, encoding='latin1')
    return(X)


def get_names(prefix):
    with gzip.open(prefix + "_attr.cp.gz", 'rb') as f:
        attr = pickle.load(f, encoding='latin1')
    return(attr['names'])

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Usage: %s *.gz" % sys.argv[0])
        sys.exit(1)
    f_list = set()
    for arg in sys.argv[1:]:
        f = arg.replace("_csr.cp.gz", '').replace("_attr.cp.gz", '')
        f_list.add(f)
    for f in f_list:
        mat = coo_matrix(load_histonemat(f, []))
        mmwrite(os.path.basename(f) + ".mtx", mat, precision=12)
