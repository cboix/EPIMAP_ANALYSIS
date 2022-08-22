#!/usr/bin/python
# ----------------------------------
# Auxiliary functions for processing
# Masterlist/other matrices
# ----------------------------------
import os
import umap
# import csv
import socket
import time
import numpy as np
import pandas as pd
# import six.moves.cPickle as pickle
from scipy.sparse import hstack, csr_matrix  # vstack

# Plotting:
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import seaborn as sns
from matplotlib import collections as mc

# Local:
import auxfunc_masterlist as AUX

rng = np.random


class process_ML(object):
    def __init__(self, mark, base, binary=True):
        self.mark = mark
        self.base = base
        self.binary = binary
        self.domain = socket.getfqdn()
        if 'broadinstitute.org' in self.domain:
            self.HOME = '/broad/compbio/cboix'
        else:
            self.HOME = os.getenv("HOME")
        self.datadir = self.HOME + '/EPIMAP_ANALYSIS/db/'
        self.mpd_dir = self.datadir + 'ChromImpute/mergedpeak_data/'
        self.types = ['imputed', 'observed']
        self.prefix = 'master_' + str(self.mark) + '_from' + str(self.base)
        self.fileprefix = self.mpd_dir + self.prefix
        self.cutoff = 2  # Need more than cutoff obs to keep a pk
        self.start = time.time()

    def setup(self):
        print("[STATUS] " + "Setting up training object")
        self.load_data()
        self.load_attr()
        self.load_info()
        self.initialize()
        self.pct_trace = []
        self.trace = []

    def get_params(self):
        params_dict = {'mark': self.mark,
                       'base': self.base}
        return(params_dict)

    def load_info(self):
        self.attr = AUX.load_pickle(self.fileprefix + "_attr.cp")
        self.pd = self.attr['peaks']
        self.cd = self.attr['cells']
        self.cells = [''] * len(self.cd)
        for k in self.cd.keys():
            self.cells[self.cd[k] - 1] = k
        self.peaks = [''] * len(self.pd)
        for k in self.pd.keys():
            self.peaks[self.pd[k] - 1] = k
        self.info = pd.read_csv(
            self.fileprefix + '_DHSs_info.tsv', sep="\t", header=None,
            names=['chr', 'start', 'end', 'chunk', 'sig',
                   'nkept', 'nkept2', 'size', 'sm1', 'sm2', 'sm3'])
        print("[STATUS] Read in " + str(self.info.shape[0]) + " infolines")

    def load_annotation(self):
        self.ann = pd.read_csv(
            self.datadir + 'Annotation/updated_imputation_metadata.tsv',
            sep="\t", header=0)
        ind = np.where([i not in ['female', 'male']
                        for i in list(self.ann.sex)])[0]
        self.ann.sex[ind] = 'unknown/mix'
        ind = np.where([i not in ['adult', 'child', 'embryonic', 'newborn']
                        for i in list(self.ann.lifestage)])[0]
        self.ann.lifestage[ind] = 'unknown/mix'
        print("[STATUS] Read in " + str(self.ann.shape[0]) + " infolines")
        self.coldf = pd.read_csv(
            self.datadir + 'Annotation/updated_tissue_colors.tsv',
            sep="\t", header=0)
        self.coldf = self.coldf.loc[list(np.argsort(self.coldf.group)), :]
        self.coldf.index = list(range(len(self.coldf)))
        self.subann = self.ann[['newgroup', 'Extended Info', 'BSSID',
                               'sex', 'lifestage', 'type']]
        self.subann.columns = ['group', 'info', 'id',
                               'sex', 'lifestage', 'type']
        self.subann = pd.merge(self.subann, self.coldf)
        self.bkdf = pd.read_csv(
            self.datadir + 'Annotation/tissue_bkpts_subsets.tsv',
            sep="\t", header=0)
        self.bk_lv = np.unique(self.bkdf.Name)
        print("[STATUS] Read in " + str(self.bkdf.shape[0]) +
              " points for dendrogram subsets")

    # Load matrix
    def load_data(self):
        self.imp = self.load_matrix('imputed')
        self.obs = self.load_matrix('observed')
        print("[STATUS] " + "Loaded data with shape (imputed): " +
              str(self.imp.shape) + " and (observed): " + str(self.obs.shape))
        self.mi = np.array(np.sum(self.imp > 0, axis=1).T)[0]
        self.mo = np.array(np.sum(self.obs > 0, axis=1).T)[0]
        self.ci = np.array(np.sum(self.imp > 0, axis=0))[0]
        self.co = np.array(np.sum(self.obs > 0, axis=0))[0]
        # Which cells are kept:
        self.ki = np.where(self.ci > 0)[0]
        self.ko = np.where(self.co > 0)[0]
        # Cutoff matrix:
        self.marg = self.mi + self.mo
        self.keptix = np.where(self.marg > self.cutoff)[0]
        self.imp = self.imp[self.keptix, ]
        self.obs = self.obs[self.keptix, ]
        print("[STATUS] " + "Keeping " + str(len(self.keptix)) +
              " peaks with evidence > " + str(self.cutoff))

    def load_matrix(self, dataset):
        filename = self.fileprefix + "_" + dataset + "_csr.cp"
        try:
            X = AUX.load_pickle(filename)
        except:
            X = AUX.load_sparse_csr(filename)
        # TODO: Fix how we load data (currently wrong indexing)
        X = X[:, 1:X.shape[1]]
        X = X[1:X.shape[0], :]
        if self.binary:
            X = 1 * (X > 0)
        return(X)

    def calc_mix_jaccard(self):
        self.X = hstack([self.imp[:, self.ki], self.obs[:, self.ko]])
        self.X = csr_matrix(self.X)
        self.pw = AUX.pairwise_jaccard(self.X.T)

    # TODO: Subset functs so we dont repeat this calculation:
    def get_mix_jaccard(self):
        if not hasattr(self, 'pw'):
            self.calc_mix_jaccard()
        self.isimp = np.array([1 for i in self.ki] + [0 for i in self.ko])
        self.collated_names = np.array([self.cells[i] for i in self.ki] +
                                       [self.cells[i] for i in self.ko])
        self.x_names = np.array(['IMP_' + self.cells[i] for i in self.ki] +
                                ['OBS_' + self.cells[i] for i in self.ko])
        # -----------------------------------------------
        # Merge info with names to get colors + infoline.
        self.colm = pd.merge(pd.DataFrame({'id': self.collated_names}),
                             self.subann, how='left')
        # Filter out non-mapping:
        self.keep = np.where([type(g) == str for g in list(self.colm.group)])[0]
        self.subX = self.X[:, self.keep]
        self.subpw = self.pw[self.keep[:, np.newaxis], self.keep]
        self.subpw[np.isnan(self.subpw)] = 0
        self.isimp = self.isimp[self.keep]
        self.x_names = self.x_names[self.keep]
        self.collated_names = self.collated_names[self.keep]
        self.colm = self.colm.loc[self.keep, :]
        # Get attributes (group # and assignments)
        self.keepcol = np.where([g in list(self.colm.group)
                                 for g in self.coldf.group])[0]
        self.dmap = {}
        self.dnum = {}
        self.d_lv = {}
        nset = 'group'
        self.d_lv[nset] = list([self.coldf.group[i] for i in self.keepcol])
        self.cc = list([self.coldf.color[i] for i in self.keepcol])
        self.dmap[nset] = plt.cm.colors.LinearSegmentedColormap.from_list(
            nset + ' colors', list(self.cc), N=len(self.cc))
        self.dnum[nset] = np.array([np.where(np.array(
            self.d_lv[nset]) == g)[0][0] for g in list(self.colm[nset])])
        # Make group and lv for categorical variables::
        self.dcol = {'sex': ['red', 'blue', 'grey'],
                     'lifestage': ['royalblue', 'purple', 'forestgreen',
                                   'indianred', 'grey'],
                     'type': ['royalblue', 'purple', 'indianred', 'goldenrod']}
        for nset in ['sex', 'lifestage', 'type']:
            self.d_lv[nset] = np.unique(list(self.ann[nset]))
            self.dmap[nset] = plt.cm.colors.LinearSegmentedColormap.from_list(
                nset + ' colors', self.dcol[nset], N=len(self.dcol[nset]))
            self.dnum[nset] = np.array([np.where(np.array(
                self.d_lv[nset]) == g)[0][0] for g in list(self.colm[nset])])

    def plot_mix_heatmap(self):
        if not hasattr(self, 'subpw'):
            print("Calculating mix")
            self.get_mix_jaccard()
        (reord, roword, colord) = AUX.reorder_by_linkage(self.subpw)
        tags_reord = [self.x_names[i] for i in roword]
        # Plot figure:
        plt.figure(figsize=(15, 15))
        sns.set(font_scale=0.25)
        ax = plt.gca()
        ax.set_facecolor('white')
        sns.heatmap(reord, cbar=False, vmax=1, vmin=0,
                    yticklabels=tags_reord, xticklabels=tags_reord)

    def make_io_lc(self):
        lines = []
        lcol = []
        for name in np.unique(self.collated_names):
            idx = np.where(self.collated_names == name)[0]
            if len(idx) == 2:
                ln = [tuple(self.u[idx[0], :]), tuple(self.u[idx[1], :])]
                col = self.cmap(self.cnum[idx[0]])
                lines.append(ln)
                lcol.append(col)
        lc = mc.LineCollection(lines, colors=lcol, linewidths=1, zorder=2)
        return(lc)

    def make_bk_lc(self):
        lines = []
        lcol = []
        for bset in self.bk_lv:
            bnames = list(self.bkdf.id[self.bkdf.Name == bset])
            # TODO: Find centroid.
            pts = []
            for name in bnames:
                idx = np.where(self.collated_names == name)[0]
                if len(idx) > 0:
                    pts.append(self.u[idx[0], :])
                    lcol.append(self.cmap(self.cnum[idx[0]]))
            pts = np.array(pts)
            ctd = np.median(pts, axis=0)
            for i in range(pts.shape[0]):
                ln = [tuple(ctd), tuple(pts[i, :])]
                lines.append(ln)
        # TODO: Filter by belongs to certain classes?
        lc = mc.LineCollection(lines, colors=lcol,
                               linestyles='dashed',
                               alpha=0.5,
                               linewidths=.8, zorder=1)
        # TODO ALSO ADD TEXT LABELS.
        return(lc)

    # MDS PLOT linking imp to observed
    def plot_mix_umap(self, metric='precomputed', show=True, size=15,
                      n_neighbors=20, min_dist=0.1, seed=None,
                      add_ioline=True, add_bkline=True, colorby='group'):
        plt.ioff()
        fit = umap.UMAP(metric=metric, random_state=seed,
                        n_neighbors=n_neighbors, min_dist=min_dist)
        if metric == 'precomputed':
            self.u = fit.fit_transform(self.subpw)
        else:
            self.u = fit.fit_transform(self.subX.T)
        # Get colors / levels:
        self.cnum = self.dnum[colorby]
        self.cmap = self.dmap[colorby]
        self.c_lv = self.d_lv[colorby]
        # Make plot:
        plt.figure(figsize=(11, 8))
        sns.set(font_scale=1)
        ax = plt.gca()
        ax.set_facecolor('white')
        iid = np.where(self.isimp == 1)[0]
        oid = np.where(self.isimp == 0)[0]
        # Concatenate the circle with an inner circle:
        circle = mpath.Path.unit_circle()
        verts = np.concatenate([circle.vertices, circle.vertices[::-1, ...]])
        codes = np.concatenate([circle.codes, circle.codes])
        open_circle = mpath.Path(verts, codes)
        ax = plt.gca()
        # Make lines collection for lines of linked observed-imputed.
        if add_ioline:
            lc = self.make_io_lc()
            ax.add_collection(lc)
        # Make lines collection for lines of subsets of breakpts:
        if add_bkline:
            lc2 = self.make_bk_lc()
            ax.add_collection(lc2)
        fig = plt.scatter(self.u[iid, 0], self.u[iid, 1], c=self.cnum[iid],
                          cmap=self.cmap, s=size, marker=open_circle, zorder=10)
        fig = plt.scatter(self.u[oid, 0], self.u[oid, 1], c=self.cnum[oid],
                          cmap=self.cmap, s=size, zorder=11)
        # title + cbar:
        plt.title('UMAP embedding of imputed (open) and observed (closed) ' +
                  self.mark)
        formatter = plt.FuncFormatter(lambda val, loc: self.c_lv[val])
        ntick = len(self.c_lv)
        cbar = plt.colorbar(ticks=list(range(ntick)), format=formatter,
                            aspect=ntick, shrink=ntick / 33)
        plt.clim(-0.5, ntick - 0.5)
        plt.tight_layout()
        cbar.ax.invert_yaxis()
        # Also plot centroids for each bkpt + label where it is.
        fig = plt.gcf()
        if show:
            plt.show()
        return(fig)

#     def write_out(self):
#         print("[STATUS] " + "Writing parameters, centers, assignments to " +
#               self.outprefix + " prefix.")
#         # Write cluster centers:
#         AUX.save_pickle(self.outprefix + "_centers.cp", self.centers)
#         # Write parameters:
#         params = self.get_params()
#         AUX.save_pickle(self.outprefix + "_params.cp", params)
#         # Plot diagnostics:
#         if self.plot_diagnostics:
#             self.plot_trace(self.outprefix + "_trace.png")
#             self.plot_centers(self.outprefix + "_centers.png")
#         # Write bedfile with assignments:
#         AUX.write_bedfile(self.outprefix + "_assignments.bed",
#                           chrom=self.keptchr, ind=self.keptix,
#                           clust=self.assign, width=200)


def hex_to_rgb(col):
    if type(col) != str:
        col = '#D3D3D3'
    if col[0] == "#" and len(col) == 7:
        col = col.lstrip('#')
        col = tuple(int(col[i:i+2], 16) for i in (0, 2, 4))
    else:
        print("Not a valid hex, returning as is")
    return(col)
