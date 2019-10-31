#!/usr/bin/python
# ----------------------------------
# Cluster large sparse matrix
# using kmeans with jaccard distance
# ----------------------------------
import os
import re
import csv
import socket
import time
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, coo_matrix  # , hstack, vstack
import fire

# Clustering:
from scipy.cluster import hierarchy as sch
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform
from numpy import linalg as LA
# from fastcluster import linkage # NOTE: BETTER.
# from polo import optimal_leaf_ordering

# Auxiliary:
import parse_config
import auxfunc_masterlist as AUX

domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')
# from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import pyplot as plt
# from IPython import display
from matplotlib import gridspec
import seaborn as sns


rng = np.random


class cluster_binary_mat(object):
    def __init__(self, filename, attrfile, outprefix, K=300, iters=None,
                 maxiters=2000, verbose=True, init_medoids=True,
                 mergestates=False, mergeset=None, extract=False,
                 h3file=None, h3attr=None, keepfile=None,
                 distmeasure='ejaccard'):
        self.filename = filename
        self.attrfile = attrfile
        self.outprefix = outprefix
        self.K = K
        self.iters = iters
        self.maxiters = maxiters
        self.verbose = verbose
        # Other arguments:
        self.init_medoids = init_medoids
        self.cutoff = 0
        self.chunksize = 100000
        # Merge:
        self.mergestates = mergestates
        self.mergeset = mergeset
        self.extract = extract
        self.keepfile = keepfile  # If keep specified indices
        # Files if merging with H3K27ac:
        self.h3file = h3file
        self.h3attr = h3attr
        self.binarymat = True
        self.distmeasure = distmeasure
        if self.distmeasure == 'euclidean':
            self.binarymat = False
        # Turn mergeset to dictionary:
        if self.mergeset is not None:
            self.mergedict = {}
            i = 0
            for ilist in self.mergeset:
                for item in ilist:
                    self.mergedict[str(item)] = str(i)
                i = i + 1
            print(self.mergedict)
        # Directories:
        self.datadir = str(os.environ['DBDIR']) + "/"
        self.keptlist = os.environ['ANNDIR'] + '/kept_bssid_20190322.txt'
        self.orderfile = os.environ['ANNDIR'] + '/bssid_order_frz20190326.tsv'
        self.chromfile = self.datadir + 'Annotation/hg19.chrom.sizes_noY'
        self.chrlist = ['chr' + str(j + 1) for j in range(22)] + ['chrX']
        self.plot_diagnostics = True
        self.closest_clust = {}
        self.pct_change = 100.0
        self.pct_cutoff = 0.01  # 0.01 in roadmap
        self.nclose = 50

    def verboseprint(self, *args):
        if self.verbose:
            for arg in args:
                print(arg)

    def main(self):
        self.setup()
        if self.extract:
            self.write_raw_regions()
        else:
            self.train()

    def setup(self):
        print("[STATUS] " + "Setting up training object")
        self.load_keptsamples()
        self.load_data()
        if self.allbin:
            self.load_chromfile()
            self.convert_ix_chrom()
        self.initialize()
        self.pct_trace = []
        self.trace = []

    def train(self):
        print("[STATUS] " + "Training...")
        self.start = time.time()
        if self.iters is None:
            for self.current_iter in range(self.maxiters):
                u, c = np.unique(self.assign, return_counts=True)
                print("COUNTS: " + str(c))
                self.step()
                if self.pct_change < self.pct_cutoff:
                    print("[STATUS] " + "Percent change below cutoff, at " +
                          str(self.pct_change) +
                          ". Writing data to " + self.outprefix)
                    self.write_out()
                    break
        else:
            for self.current_iter in range(self.iters):
                self.step()
            self.write_out()

    def step(self):
        self.t2assign = self.oldassign.copy()
        self.oldassign = self.assign.copy()
        self.calc_distance()  # And re-assign
        self.calc_centers()
        # Timing:
        self.current_time = time.time() - self.start
        print
        print("[STATUS] " + str(self.current_iter) + ": " +
              str(round(self.current_time, 2)) + "s")
        # Print and save number changed
        self.changes = (self.oldassign - self.assign) != 0
        self.changeloc = np.where(self.changes)[0]
        print("Changes in: " + str(self.changeloc[0:10]))
        self.pct_change = (np.sum(self.changes) * 1.0 / self.M) * 100
        print(list(self.oldassign[0:10]))
        print(list(self.assign[0:10]))
        # Count changes between old + new in last 10 iter:
        # Locations that reverted:
        self.fliploc = [i for i in self.changeloc
                        if self.t2assign[i] == self.assign[i]]
        self.flips[self.fliploc] = self.flips[self.fliploc] + 1
        print("Number of flips: " + str(len(self.fliploc)))
        self.pct_trace.append(self.pct_change)
        self.trace.append(self.cost)
        print("[STATUS] " + str(round(self.pct_change, 3)) +
              "%" + " assignments changed.")
        if self.current_iter % 10 == 0:
            self.write_out()

    def get_params(self):
        params_dict = {
            'filename': self.filename,
            'K': self.K, 'pct_trace': self.pct_trace,
            'trace': self.trace,
            'last_iter': self.current_iter,
            'cutoff': self.cutoff,
            'time': self.current_time}
        return(params_dict)

    def write_out(self):
        print("[STATUS] " + "Writing parameters, centers, assignments to " +
              self.outprefix + " prefix.")
        # Write cluster centers:
        print("Centers...")
        AUX.save_pickle(self.outprefix + "_centers.cp", self.centers)
        np.savetxt(self.outprefix + "_centers.tsv.gz",
                   self.centers, delimiter="\t", fmt='%1.4f')
        # Write parameters:
        print("Params...")
        params = self.get_params()
        AUX.save_pickle(self.outprefix + "_params.cp", params)
        # Plot diagnostics:
        if self.plot_diagnostics:
            print("Plots...")
            self.plot_trace(self.outprefix + "_trace.png")
            self.plot_centers(self.outprefix + "_centers.png")
            # self.plot_matrix(self.outprefix + "_matrices.png")
        # Write bedfile with assignments:
        self.outfile = self.outprefix + "_assignments.bed"
        print("Writing assignments to: " + self.outfile)
        if self.allbin:
            write_bedfile(self.outfile, chrom=self.keptchr, ind=self.keptix,
                          clust=self.assign, width=200)
        else:
            # Use pandas to_csv to write faster:
            # NOTE: Write with 1-indexing:
            df = pd.DataFrame({'id': self.keptix + 1,
                               'cls': self.assign})
            df.to_csv(self.outfile, sep='\t', index=False)

    def write_raw_regions(self):
        print("[STATUS] Writing " + str(self.X.nnz) +
              " nnz elements in X matrix.")
        Xcoo = self.X.tocoo()
        # Outprefix serves as directory as well:
        print(self.X.shape[1])
        for i in range(self.X.shape[1]):
            self.outfile = self.outprefix + "/raw_%03d_assignments.bed" % i
            ind = Xcoo.row[Xcoo.col == i]
            outind = self.keptix[ind]
            print(i, len(outind))
            df = pd.DataFrame({'id': outind + 1, 'cls': i})
            df.to_csv(self.outfile, sep='\t', index=False)

    # Load matrix
    def load_data(self):
        try:
            print("[STATUS] Trying to load cp gz: " + self.filename)
            self.X = AUX.load_pickle_gzip(self.filename)
        except:
            print("[STATUS] Trying to load sparse csr: " + self.filename)
            self.X = load_sparse_csr(self.filename)
        # Load attributes:
        print("[STATUS] Load attributes in cp gz: " + self.attrfile)
        self.attr = AUX.load_pickle_gzip(self.attrfile)
        self.names = self.attr['names']
        self.names_split = [n.split("_")[0] for n in self.names]
        try:
            self.states_split = [n.split("_")[1] for n in self.names]
        except:
            print("[STATUS] Do not need to merge states because already merged")
            self.states_split = []
            self.mergestates = False
        # Collapse matrix if necessary:
        if self.mergestates:
            self.collapse_matrix()
        # Merge after collapse/not:
        if self.h3file is not None:
            self.merge_matrices(binary=self.binarymat)
        # Get plotting order:
        self.nidx = []
        for n in self.names_split:
            self.nidx.append([i for i, v in enumerate(self.order) if v == n][0])
        print("Saving rownames...")
        AUX.write_list(self.outprefix + "_names.tsv", self.names)
        print("NNAMES: " + str(len(self.names)))
        # Collapse data:
        self.allbin = (self.X.shape[0] > 15000000)
        print("[STATUS] " + "Loaded data with shape: " + str(self.X.shape))
        if self.keepfile is None:
            self.marg = np.sum(self.X, axis=1)
            # Locate points above cutoff
            self.keptix = np.where(self.marg > self.cutoff)[0]
            # NOTE: For future versions, merge peaks, not straighforward now
            # Keep locations for original peaks
            print("[STATUS] " + "Keeping " + str(len(self.keptix)) +
                  " peaks with evidence > " + str(self.cutoff))
        else:
            print("[STATUS] Loading indices from: " + self.keepfile)
            self.keptix = AUX.load_pickle_gzip(self.keepfile)
            print("[STATUS] Keeping " + str(len(self.keptix)) + " peaks")
        self.X = self.X[self.keptix, ] * 1.0
        self.M = self.X.shape[0]
        self.N = self.X.shape[1]

    def collapse_matrix(self):
        print("[STATUS] Collapsing matrix")
        # 1. Make mapping from old to new indices
        namedict = {}
        j = 0
        fromidx = []
        toidx = []
        collapsed_names = []
        if self.mergeset is None:
            # Merge all with same name
            for i in range(len(self.names)):
                nam = self.names_split[i]
                if nam not in namedict.keys():
                    namedict[nam] = j
                    collapsed_names.append(nam)
                    j = j + 1
                fromidx.append(i)
                toidx.append(namedict[nam])
        else:
            # Merge according to mergeset groups
            for i in range(len(self.names)):
                nam = self.names_split[i]
                st = self.states_split[i]
                st_tag = self.mergedict[st]
                fulltag = nam + "_" + st_tag
                if fulltag not in namedict.keys():
                    namedict[fulltag] = j
                    collapsed_names.append(fulltag)
                    j = j + 1
                fromidx.append(i)
                toidx.append(namedict[fulltag])
        # 2. Collapse the matrix
        Xcoo = self.X.tocoo()
        toidx = np.array(toidx)
        newcol = toidx[Xcoo.col]
        # 3. Replace X matrix with collapsed mat:
        newcoo = coo_matrix((Xcoo.data, (Xcoo.row, newcol)),
                            shape=(Xcoo.shape[0], j))
        self.X = newcoo.tocsr()
        print("[STATUS] Collapsed matrix to: " + str(self.X.shape))
        # Update names and names_split:
        self.names = collapsed_names
        self.names_split = [n.split("_")[0] for n in self.names]

    def merge_matrices(self, binary=False):
        try:
            print("[STATUS] Trying to load cp gz: " + self.h3file)
            self.H = AUX.load_pickle_gzip(self.h3file)
        except:
            print("[STATUS] Trying to load sparse csr: " + self.h3file)
            # NOTE: Won't work unless strip off .cp.gz:
            self.H = load_sparse_csr(self.h3file)
        # Load attributes:
        print("[STATUS] Load attributes in cp gz: " + self.h3attr)
        self.hattr = AUX.load_pickle_gzip(self.h3attr)
        self.hnames = self.hattr['names']
        tmpnames = [re.sub("^.*=BSS", "BSS", n) for n in self.hnames]
        self.hnames = [re.sub("_.*$", "", n).strip("\n") for n in tmpnames]
        print(self.hnames[0:10])
        # NOTE: Set 0 if lower than cutoff=2
        cutoff = 2.0
        if self.X.shape == self.H.shape:
            # Make sure that rows etc match:
            hord = np.array([np.where(np.array(self.hnames) == n)[0][0]
                             for n in self.names])
            print("[STATUS] Reordering histone matrix")
            print(self.names[0:10])
            print(self.hnames[0:10])
            print(hord[0:10])
            # Turn hX into binary matrix, thresholded at 2:
            self.Hord = 1 * (self.H[:, hord] >= cutoff)
            del(self.H)
            Xtmp = self.X.multiply(self.Hord)
            print(str(self.X.nnz) + " (ann) to " + str(Xtmp.nnz))
            print(str(self.Hord.nnz) + " (histone) to " + str(Xtmp.nnz))
            self.X = Xtmp
            del(Xtmp)
        else:
            Hcoo = self.H.tocoo()
            print("[STATUS] Reducing to values above " + str(cutoff))
            ind = np.where(Hcoo.data >= cutoff)[0]
            Hred = coo_matrix((Hcoo.data[ind], (Hcoo.row[ind], Hcoo.col[ind])),
                              shape=self.H.shape)
            print(str(Hcoo.nnz) + " to " + str(Hred.nnz))
            del(self.H, Hcoo)
            self.H = Hred.tocsc()
            del(Hred)
            # Arrays for new coo matrix:
            Xcsc = self.X.tocsc()
            for i in range(len(self.hnames)):
                nam = self.hnames[i]
                print(i, nam)
                # If name in list:
                mtid = [j for j, v in enumerate(self.names_split)
                        if v == nam]
                print(mtid)
                if len(mtid) > 0:
                    mtid = np.array(mtid)
                    vec = self.H[:, i]  # Turn to vec
                    print(vec.shape)
                Xcsc[:, mtid] = Xcsc[:, mtid].multiply(vec)
            print(str(self.X.nnz) + " to " + str(Xcsc.nnz))
            self.X = Xcsc.tocsr()
            del(Xcsc, self.H)

    def load_keptsamples(self):
        self.samples = []
        self.order = []
        with open(self.keptlist, 'r') as f:
            for line in f:
                self.samples.append(line.split("\n")[0])
        with open(self.orderfile, 'r') as f:
            for line in f:
                samp = line.split("\t")[0]
                if samp in self.samples:
                    self.order.append(samp)

    def load_chromfile(self):
        with open(self.chromfile, 'r') as fp:
            reader = csv.reader(fp, delimiter='\t')
            self.chromsizes = {}
            for line in reader:
                self.chromsizes[line[0]] = int(line[1])

    # Find corresponding chrom for indexes:
    def convert_ix_chrom(self):
        print("[STATUS] " + "Finding chromosomes corresponding " +
              "to matrix indices.")
        last = 0
        self.chr_idx = []
        for chrom in self.chrlist:  # Iter in order.
            size = int(np.floor(self.chromsizes[chrom] / 200.0))
            last = last + size
            self.chr_idx = self.chr_idx + [chrom] * size
        self.chr_idx = np.array(self.chr_idx)
        self.keptchr = self.chr_idx[self.keptix]
        print("[STATUS] Read in and created " + str(last) + " indices. " +
              "Length of chr idx vector: " + str(len(self.chr_idx)))

    # Initialize variables:
    def initialize(self):
        if self.init_medoids:
            # Init clusters with random points:
            self.cid = rng.choice(self.M, self.K, replace=False)
            self.centers = self.X[self.cid, :].toarray().T * 1.0
        else:
            # Init clusters with random assignments using xavier init:
            inner = rng.randn(self.N, self.K) * 2
            self.centers = 1.0 / (1.0 + np.exp(-inner))
        self.center_marg = np.squeeze(np.array(np.sum(self.centers, axis=0)))
        print(self.center_marg)
        # Margin for distance calculations:
        self.X_marg = np.squeeze(np.array(np.sum(self.X, axis=1)))
        # Space allocation for assignment:
        self.assign = np.zeros(self.M).astype(int)
        self.oldassign = self.assign.copy()
        self.flips = self.assign.copy()

    # Calculate cluster centers from weights
    def calc_centers(self):
        for i in range(self.K):
            ix = np.where(self.assign == i)[0]
            if (len(ix) > 0):
                out = np.squeeze(np.array(np.sum(self.X[ix, ], axis=0) /
                                          (1.0 * len(ix))))
                self.centers[:, i] = out
        self.center_marg = np.squeeze(np.array(np.sum(self.centers, axis=0)))

    # Get closest from distance between centers:
    def get_closest_clusters(self):
        self.c_dist = squareform(pdist(self.centers.T, metric='cosine'))
        for i in range(self.K):
            self.closest_clust[i] = np.argsort(self.c_dist[:, i])[0:self.nclose]

    # NOTE: May be faster to make assign a sparse matrix:
    # Of size M x K and multiply centers = X.T.dot(assign),
    # followed by normalizing norm by size.
    def calc_distance(self):
        self.get_closest_clusters()
        if self.pct_change < 10.0:
            self.verboseprint("[STATUS] Calculating reduced distance")
            self.red_calc_distance()
        else:
            self.verboseprint("[STATUS] Calculating full distance")
            self.calc_distance_by_chunks()

    # Calculate distance from all cluster centers
    def calc_distance_by_chunks(self):
        self.cost = 0
        NCHUNKS = int(self.M / self.chunksize) + 1
        for i in range(NCHUNKS):
            ix = range((i * self.chunksize),
                       np.min([self.M, (i+1) * self.chunksize]))
            # print(ix)
            # fJaccard: ( sum_min_xy / sum_max_xy )
            if self.distmeasure == 'euclidean':
                dist = LA.norm(self.X[ix, ], self.centers)
            else:
                dist = jaccard_dist(self.X[ix, ], self.centers,
                                    self.X_marg[ix, None],
                                    self.center_marg,
                                    jitter=True)
            self.assign[ix] = np.argmin(dist, 1)
            self.cost = self.cost + np.sum(np.min(dist, 1))

    # Option 1: Look at closest
    # Option 2: Do Minibatch
    # Option 3: Do both.
    def red_calc_distance(self):
        # Keep tmp, changing assignments:
        tmpassign = self.assign.copy()
        for i in range(self.K):
            # Pick corresponding to cls i
            ix = np.where(self.assign == i)[0]
            if len(ix) > 0:
                closest = self.closest_clust[i]
                reduced_center = self.centers[:, closest]
                if self.distmeasure == 'euclidean':
                    dist = LA.norm(self.X[ix, ], reduced_center)
                else:
                    dist = jaccard_dist(self.X[ix, ], reduced_center,
                                        self.X_marg[ix, None],
                                        self.center_marg[closest],
                                        jitter=True)
                # Assign min, mapping through closest set:
                tmpassign[ix] = closest[np.argmin(dist, 1)]
        self.assign = tmpassign

    def plot_trace(self, filename):
        sns.set(font_scale=1.1)
        # Plot current trace:
        fig = plt.figure(figsize=(9, 12))
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
        ax = plt.subplot(gs[0])
        # Plot % change:
        ax.set_title("Percent assignments changing")
        ax.set_facecolor('white')
        y = np.array(self.pct_trace)
        if len(y) > 0:
            x = range(len(y))
            plt.plot(x, y, lw=2)
        # Plot trace:
        ax = plt.subplot(gs[1])
        ax.set_title("Trace: Cost")
        ax.set_facecolor('white')
        y = np.array(self.trace[1:])
        if len(y) > 0:
            x = range(len(y))
            plt.plot(x, y, lw=2)
        plt.tight_layout()
        fig = plt.gcf()
        fig.savefig(filename, dpi=150, bbox_inches='tight')

    def plot_centers(self, filename, diagonal=True, reordsamp=False):
        # Counts (+ reordered):
        u, c = np.unique(self.assign, return_counts=True)
        c = c[np.argsort(u)]
        # Reorder matrix, using only vals that exist:
        # centers = N x K (cells by types). NOTE: cols = kept clusters
        if reordsamp:
            mat = self.centers[:, np.sort(u)].T
            [reord, row_idx, col_idx] = reorder_by_linkage(mat)
            if diagonal:
                (reord, factor_ord, cols) = diagonalize_mat(
                    mat[:, col_idx], cutoff=0.05)
                # Cols tells us which to ignore:
                c = c[cols]
                c = c[factor_ord]
        else:
            # Reorder mat by names, using fixed order (to idx, then order):
            self.nameord = np.argsort(self.nidx)
            mat = self.centers[self.nameord, :]
            mat = mat[:, np.sort(u)].T
            if diagonal:
                (reord, factor_ord, cols) = diagonalize_mat(mat, cutoff=0.05)
                # Cols tells us which to ignore:
                c = c[cols]
                c = c[factor_ord]
        col_lines = np.arange(5, reord.shape[1], 5)
        # -------------------
        sns.set(font_scale=1.5)
        plt.figure(figsize=(8, 15))
        gs = gridspec.GridSpec(2, 2, height_ratios=[2, .5],
                               wspace=0.01, hspace=0.01)
        ax = plt.subplot(gs[0])
        sns.heatmap(reord.T, cmap='Blues', ax=ax, cbar=False,
                    xticklabels=False, yticklabels=False, vmin=0, vmax=1)
        plt.ylabel('Cell Types/States')
        ax.vlines(col_lines, *ax.get_ylim(), lw=.5, alpha=.5,
                  linestyle='dashed', color='grey')
        # ---------------
        ax = plt.subplot(gs[1])
        sns.heatmap(1 * (reord.T > 0.1), cmap='Blues', ax=ax, cbar=False,
                    xticklabels=False, yticklabels=False, vmin=0, vmax=1)
        ax.vlines(col_lines, *ax.get_ylim(), lw=.5, alpha=.5,
                  linestyle='dashed', color='grey')
        # ---------------
        ax = plt.subplot(gs[2])
        ind = range(len(c))
        ax.set_xlim([- 0.5, len(c) - .5])
        ax.set_ylim([0, np.max(c)])
        ax.set_facecolor('white')
        ax.bar(ind, c)
        plt.xlabel('Clusters')
        ax.vlines(col_lines - 0.5, *ax.get_ylim(), lw=.5, alpha=.5,
                  linestyle='dashed', color='grey')
        ax.invert_yaxis()
        # Save fig:
        plt.tight_layout()
        fig = plt.gcf()
        fig.savefig(filename, dpi=450, bbox_inches='tight')

    def plot_matrix(self, filename, diagonal=True):
        # Plot the matrix (X)
        # Using imshow or spy
        sns.set(font_scale=1.5)
        plt.figure(figsize=(8, 20))
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace=0.05)
        ax = plt.subplot(gs[0])
        ax.set_facecolor('white')
        # ax.spy(self.X)  # cmap=plt.cm.binary)
        tiny_spy(self.X, ax)
        # self.Xdense = self.X.toarray()
        # plt.imshow(self.Xdense, cmap=plt.cm.binary)
        # Plot the matrix reordered according to assignments.
        ax = plt.subplot(gs[1])
        ax.set_facecolor('white')
        # neword = np.argsort(self.assign)
        # ax.spy(self.X[neword, :])  # cmap=plt.cm.binary)
        # tiny_spy(self.X, ax, rowreord=neword)
        plt.ylabel('Cell Types/States')
        plt.xlabel('Locations')
        # Save fig:
        fig = plt.gcf()
        fig.savefig(filename, dpi=350)


def diagonalize_mat(mat, cutoff=0.05):
    Z = mat > cutoff
    colsums = np.sum(Z, axis=1)
    cols = np.array([i for i, v in enumerate(colsums) if v > 0])
    maxids = [np.max(np.where(Z[i, :])[0]) for i in cols]
    diag_idx = np.argsort(maxids)
    reord = mat[cols, :]
    reord = reord[diag_idx, ]
    return(reord, diag_idx, cols)


def pairwise_jaccard(X):
    # X = X.astype(bool).astype(int)
    intrsct = X.dot(X.T)
    # row_sums = np.squeeze(np.array(np.sum(X, axis=1)))
    row_sums = intrsct.diagonal()
    unions = row_sums[:, None] + row_sums - intrsct
    dist = 1.0 - intrsct / unions
    return(dist)


def jaccard_dist(A, B, Am, Bm, jitter=False):
    intersect = A.dot(B)
    unions = Am + Bm - intersect
    if jitter:
        intjitter = intersect + rng.random(intersect.shape) * 1e-8
        dist = 1.0 - intjitter / unions
    else:
        dist = 1.0 - intersect / unions
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


# Note procedure to speed up would be to consider only
# the closest couple of centers. -> cluster centers, and choose only centers
# close to it.


def write_bedfile(filename, chrom, ind, clust, width=200):
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


def load_sparse_csr(filename):
    loader = np.load(filename + '.npz')
    return csr_matrix((loader['data'], loader['indices'], loader['indptr']),
                      shape=loader['shape'])


def custom_spy(m, ax, rowreord=None):
    from scipy.sparse import coo_matrix
    from matplotlib.patches import Rectangle
    if not isinstance(m, coo_matrix):
        m = coo_matrix(m)
    if rowreord is not None:
        print("reordering rows")
        ylist = [rowreord[i] for i in m.row]
    else:
        ylist = m.row
    for (x, y) in zip(m.col, ylist):
        ax.add_artist(Rectangle(
            xy=(x-0.5, y-0.5), width=1, height=1))
    ax.set_xlim(-0.5, m.shape[1]-0.5)
    ax.set_ylim(-0.5, m.shape[0]-0.5)
    ax.invert_yaxis()
    ax.set_aspect(float(m.shape[0])/float(m.shape[1]))


def tiny_spy(m, ax, rowreord=None):
    from scipy.sparse import coo_matrix
    if not isinstance(m, coo_matrix):
        m = coo_matrix(m)
    if rowreord is not None:
        print("reordering rows")
        ylist = [rowreord[i] for i in m.row]
    else:
        ylist = m.row
    print("plotting")
    plt.scatter(m.col, ylist, marker='.')
    ax.set_xlim(-0.5, m.shape[1]-0.5)
    ax.set_ylim(-0.5, m.shape[0]-0.5)
    ax.invert_yaxis()
    # ax.set_aspect(float(m.shape[0])/float(m.shape[1]))


if __name__ == "__main__":
    fire.Fire(cluster_binary_mat)
