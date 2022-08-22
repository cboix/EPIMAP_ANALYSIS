#!/usr/bin/python
# -------------------------------------------
# Plot a large sparse matrix using datashader
# NOTE: Plot as points first
# TODO: how to plot as properly sized pts if zoom in
# -------------------------------------------
import os
import fire
# import h5py
import time
import six.moves.cPickle as pickle
from scipy.sparse import coo_matrix, csr_matrix
import gzip

import parse_config

# For plotting points
import datashader.spatial.points as dsp
import pandas as pd
import numpy as np
import datashader as ds
import datashader.transfer_functions as tf
from collections import OrderedDict as odict

from matplotlib import pyplot as plt
from matplotlib import gridspec
import seaborn as sns
# Clustering:
from scipy.cluster import hierarchy as sch
from scipy.cluster.hierarchy import linkage
from polo import optimal_leaf_ordering

# For interactive:
import holoviews as hv
from holoviews.operation.datashader import datashade
hv.extension('bokeh')

# For random data:
rng = np.random


# Load X from file or give obj
class ds_mat(object):
    def __init__(self, filepref, reord=True, pqt=True, filetype='binary',
                 clsord=None, frozen_order=True, sampord=None):
        self.filepref = filepref
        self.reord = reord
        self.pqt = pqt
        self.frozen_order = False
        self.sampord = sampord
        self.filetype = filetype
        self.datafile = filepref + "_csr.cp.gz"
        self.attrfile = filepref + "_attr.cp.gz"
        self.clsord = clsord
        # Directories:
        self.datadir = os.environ['HOME'] + os.environ['DBDIR']
        print(self.datadir)
        self.pqloc = filepref + '_sorted.parq'
        self.assignfile = self.datadir + '/cls_300_assignments.bed'
        self.keptlist = self.datadir + '/Annotation/kept_bssid_20190322.txt'
        self.orderfile = self.datadir + '/Annotation' + \
            '/bssid_order_frz20190326_samponly.tsv'
        self.start = time.time()

    def main(self):
        self.process()
        self.plot_bycls()

    def process(self):
        self.load_data()
        self.load_metadata()
        if self.reord:
            self.reorder_data()
        self.to_dataframe()
        if self.pqt:
            if not os.path.exists(self.pqloc):
                self.to_parquet()
            self.read_parquet()

    def load_data(self):
        # Read in matrix:
        with gzip.open(self.datafile, 'rb') as f:
            self.X = pickle.load(f, encoding='latin1')

    def load_metadata(self):
        # Kept samples and frozen order:
        self.samples = []
        self.order = []
        with open(self.keptlist, 'r') as f:
            for line in f:
                self.samples.append(line.split("\n")[0])
        # Order:
        with open(self.orderfile, 'r') as f:
            for line in f:
                samp = line.split("\n")[0]
                if samp in self.samples:
                    self.order.append(samp)
        # Names:
        with gzip.open(self.attrfile, 'rb') as f:
            outattr = pickle.load(f, encoding='latin1')
        self.names = outattr['names']
        if self.filetype == 'bw':
            nout = [n.split("name=")[1] for n in self.names]
            self.names_split = [n.split("_")[0] for n in nout]
        else:
            self.names_split = [n.split("_")[0] for n in self.names]
        # For reordering:
        self.nidx = []
        for n in self.names_split:
            self.nidx.append([i for i, v in enumerate(self.order)
                              if v == n][0])
        # Read in assignments:
        self.xloc = []
        self.assign = []
        with open(self.assignfile, 'r') as f:
            for line in f:
                line = line.strip("\n").split("\t")
                self.xloc.append(int(line[0]))
                self.assign.append(int(line[1]))
        self.assign = np.array(self.assign)
        self.xloc = np.array(self.xloc)

    def reorder_cols(self, mat, method='ward'):
        print("Reordering columns by jaccard distance")
        tM = csr_matrix(mat.T) * 1.0
        [self.dist, self.intrsct] = pairwise_jaccard(tM)
        Z = linkage(self.dist, method=method)
        Z = optimal_leaf_ordering(Z, self.dist)
        col_order = sch.leaves_list(Z)
        return(col_order)

    def reorder_data(self):
        # Subset used rows and reorder cols too according to frozen order:
        if self.frozen_order:
            self.sampord = np.argsort(self.nidx)
        else:
            if self.sampord is None:
                self.sampord = self.reorder_cols(self.X)
        Xtmp = self.X[self.xloc, :]
        self.Xsord = Xtmp[:, self.sampord]
        del(Xtmp)
        print(str(self.X.shape) + " to " + str(self.Xsord.shape))
        # Calculate and diagonalize centers matrix:
        ua, c = np.unique(self.assign, return_counts=True)
        self.centers = np.zeros((self.X.shape[1], len(ua)))
        for i in ua:
            slc = np.where(self.assign == i)[0]
            self.centers[:, i] = np.array(np.mean(self.Xsord[slc, :], 0))[0]
        Zcent = 1 * (self.centers > 0.1)
        self.bymean = np.argsort(np.array([np.min(np.where(x > 0)[0])
                                           for x in Zcent.T]))
        # self.bymean = np.argsort(np.array([np.mean(np.where(x > 0)[0])
        #                                    for x in Zcent.T]))
        # Reorder assigned clusters according to clsord:
        if self.clsord is None:
            self.clsord = self.bymean
        self.meanmap = np.argsort(self.clsord)
        ordassign = [self.meanmap[i] for i in self.assign]
        self.aind = np.argsort(ordassign)
        self.sortedassign = self.assign[self.aind]

    def plot_centers(self):
        # Confirm centers are ok:
        fig = plt.figure(figsize=(10, 10))
        fig = sns.heatmap(np.matrix(self.centers[:, self.clsord]) > 0.1,
                          cmap=plt.cm.binary)
        return(fig)

    def plot_margins(self, horiz=True):
        occe = np.array(np.sum(self.Xsord, axis=1).T)[0]
        nume = np.array(np.sum(self.Xsord, axis=0))[0]
        if horiz:
            fig = plt.figure(figsize=(15, 5))
            gs = gridspec.GridSpec(1, 2)
        else:
            fig = plt.figure(figsize=(8, 10))
            gs = gridspec.GridSpec(2, 1)
        ax = plt.subplot(gs[0])
        h = plt.hist(occe, 50, color='darkgrey')
        plt.xlabel('Number of epigenomes with enhancer')
        plt.ylabel('Number of enhancer')
        ax = plt.subplot(gs[1])
        h = plt.hist(nume, 30, color='darkgrey')
        del(h)
        plt.xlabel('Number of enhancers in epigenome')
        plt.ylabel('Number of epigenomes')
        return(fig)

    def to_dataframe(self):
        print("Turning sparse matrix to plotting dataframe")
        # Subset to kept:
        if self.reord:
            self.Xred = self.Xsord[self.aind, :]
            self.Xcoo = coo_matrix(self.Xred)
        else:
            self.Xcoo = coo_matrix(self.X)
        print(self.Xcoo.shape)
        # Cutoff:
        self.xlim = (0, self.Xcoo.shape[0])
        self.ylim = (0, self.Xcoo.shape[1])
        if self.filetype == 'bw':
            # Filter:
            # TODO: CHANGE TO >= when no longer integer data!
            self.cind = self.Xcoo.data > 2
            tmpcoo = coo_matrix((self.Xcoo.data[self.cind],
                                 (self.Xcoo.row[self.cind],
                                 self.Xcoo.col[self.cind])),
                                shape=self.Xcoo.shape)
            self.Xcoo = tmpcoo
        if self.reord:
            self.srtamt = self.sortedassign[self.Xcoo.row]
        # To dataframe:
        self.xdf = pd.DataFrame(odict([('x', self.Xcoo.col),
                                       ('y', self.Xcoo.row),
                                       ('v', self.Xcoo.data)]))
        if self.reord:
            self.xclsdf = pd.DataFrame(odict([('x', self.Xcoo.col),
                                              ('y', self.Xcoo.row),
                                              ('v', self.Xcoo.data),
                                              ('cls', self.srtamt)]))
            self.xclsdf["cls"] = self.xclsdf["cls"].astype("category")
        print(self.xdf.shape)

    def to_parquet(self):
        startpq = time.time()
        dsp.to_parquet(self.xdf, self.pqloc, 'x', 'y',
                       shuffle='disk', npartitions=32)
        print("Time to parquet: " + str(round(time.time() - startpq, 2)) + "s")

    def read_parquet(self):
        # Read as spatial (subset operations are faster):
        self.sframe = dsp.read_parquet(self.pqloc)
        print(self.sframe)
        print('x_range: ', self.sframe.spatial.x_range)
        print('y_range: ', self.sframe.spatial.y_range)

    def plot_matrix(self, horiz=True, thresh=None):
        if hasattr(self, 'sframe'):
            df = self.sframe.spatial_query(self.ylim, self.xlim)
        else:
            df = self.xdf
        # Threshold:
        if thresh is not None:
            odf = df[df.v > thresh]
            fig = tf.Image(create_image(odf, self.xlim, self.ylim, horiz=horiz))
        else:
            fig = tf.Image(create_image(df, self.xlim, self.ylim, horiz=horiz))
        return(fig)

    def plot_percentiles(self):
        h = 600
        canvas = ds.Canvas(x_range=self.xlim,
                           y_range=self.ylim,
                           plot_width=h * 5, plot_height=h)
        agg = canvas.points(self.xdf, 'y', 'x')
        # fig = tf.Images([tf.shade(agg.where(agg >= np.percentile(agg, pct)),
        #                           name=str(pct) + "th Percentile")
        #                  for pct in [25, 50, 75, 90, 98, 99]])
        fig = tf.Images(
            tf.shade(agg.where(agg >= np.percentile(agg, 25)),
                     name="25th Percentile"),
            tf.shade(agg.where(agg >= np.percentile(agg, 50)),
                     name="50th Percentile"),
            tf.shade(agg.where(agg >= np.percentile(agg, 75)),
                     name="75th Percentile"),
            tf.shade(agg.where(agg >= np.percentile(agg, 90)),
                     name="90th Percentile"),
            tf.shade(agg.where(agg >= np.percentile(agg, 98)),
                     name="98th Percentile"),
            tf.shade(agg.where(agg >= np.percentile(agg, 99)),
                     name="99th Percentile"))
        fig.num_cols = 1
        # fig = plt.gcf()
        return(fig)

    def plot_bycls(self, thresh=None):
        h = 600
        pal = ds.colors.Set1
        pal = ['darkblue', 'firebrick']
        emap = pal * int(len(self.meanmap) / len(pal) + 1)
        self.extcmap = [emap[i] for i in self.meanmap]
        canvas = ds.Canvas(x_range=self.xlim,
                           y_range=self.ylim,
                           plot_width=h * 5, plot_height=h)
        if thresh is not None:
            odf = self.xclsdf[self.xclsdf.v > thresh]
            agg = canvas.points(odf, 'y', 'x', ds.count_cat('cls'))
        else:
            agg = canvas.points(self.xclsdf, 'y', 'x', ds.count_cat('cls'))
        fig = tf.Image(tf.shade(agg, color_key=self.extcmap),
                       name="With clusters")
        return(fig)

    def plot_interactive(self):
        fig = datashade(hv.Points(self.sframe, kdims=['y', 'x']),
                        cmap=plt.cm.binary,
                        min_alpha=0).opts(width=2400, height=600)
        return(fig)


# NOTE: May speed up by giving csr_matrix
def pairwise_jaccard(X):
    intrsct = X.dot(X.T)
    # row_sums = np.squeeze(np.array(np.sum(X, axis=1)))
    row_sums = intrsct.diagonal()
    unions = row_sums[:, None] + row_sums - intrsct
    dist = 1.0 - intrsct / unions
    return(dist, intrsct)


def create_image(df, x_range, y_range, horiz=True):
    scale = 5
    h = 600
    w = h * scale
    if horiz:
        canvas = ds.Canvas(x_range=x_range, y_range=y_range,
                           plot_width=w, plot_height=h)
        agg = canvas.points(df, 'y', 'x')
    else:
        canvas = ds.Canvas(x_range=y_range, y_range=x_range,
                           plot_width=h, plot_height=w)
        agg = canvas.points(df, 'x', 'y')
    return(tf.shade(agg))


if __name__ == "__main__":
    fire.Fire(ds_mat)
