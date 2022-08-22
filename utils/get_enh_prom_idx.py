#!/usr/bin/python
# ---------------------------------------------------
# Get the enhancer and promoter indices for a dataset
# ---------------------------------------------------
# Run as:
# conda activate pytorch_env > /dev/null; # has tqdm
# python ${BINDIR}/get_enh_prom_idx.py main --mixobs --nonovl --makeplots --usemark H3K27ac --savematrices
# python ${BINDIR}/get_enh_prom_idx.py main --mixobs --nonovl --makeplots --usemark None --savematrices
# python ${BINDIR}/get_enh_prom_idx.py main --mixobs --makeplots --usemark H3K27ac --savematrices
# python ${BINDIR}/get_enh_prom_idx.py main --mixobs --makeplots --usemark None --savematrices

# System:
import os
import re
import csv
import gzip
import sys
import gc
import fire
import six.moves.cPickle as pickle
from tqdm import tqdm
from scipy.io import mmwrite

# Math
import numpy as np
import pandas as pd

# Plotting:
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib import gridspec, ticker, markers

# Local paths and functions:
import auxfunc_masterlist as AUX
import parse_config

# For random data:
rng = np.random

# Seaborn style
# sns.set(style=ticks, color_codes=True)

# Statistics of ML intersections
#     #### Statistics on occupancy, intersection with the masterlist:
#     1. Number of locations labeled Enhancer, Promoter, or Both
#     2. Distributions
#         - How many enh/prom per epigenome (XXX +/- YYY)
#         - How many epigenome per enh/prom
#     3. Saturation curve (rarefaction)
#         - Starting with the matched roadmap samples
#         - Overall
#         - Which epigenomes have the larges # of "uncommon elements" - aka - potentially need more study
#     4. How much redundancy in the masterlist?

### Defining promoters/enhancers/dyadic:
# - **Roadmap approach:** We arbitrate on these regions by first clustering them (using the methods in the following section) with an expected cluster size of 10,000 regions, and then for each cluster calculating (a) the mean posterior probabilities for promoter and enhancer calls separately, and (b) the mean number of reference epigenomes in which regions were called promoter or enhancer. Clusters of regions for which the differences in mean posterior probabilities (a) is smaller than 0.05, or for which the absolute log 2 ratio of the number of epigenomes called as promoter or enhancer (b) is smaller than 0.05, are called true ‘dyadic’ regions, along with a small number of ‘ambiguous’ regions in state BivFlnk. Note that this particular clustering is only to arbitrate on these regions using group statistics instead of one-by-one;
# - **Approach here** Based on statistics, choose elements as regions with 75%+ of that element. Ignores 122k regions -> dyadic


class calculate_functional_idx(object):
    def __init__(self, mixobs=True, nonovl=False, model=None, usemark='H3K27ac',
                 makeplots=False,
                 savematrices=False, saveraw=False, savemtx=False):
        self.mixobs = mixobs
        self.nonovl = nonovl
        self.usemark = usemark
        self.makeplots = makeplots
        self.savematrices = savematrices
        self.saveraw = saveraw
        self.savemtx = savemtx
        # Set up data locations:
        self.statesuf = '_r200_e0_allchr_merged'
        self.marksuf = '_r25_e100_allchr_merged'
        if self.nonovl:
            self.datasuf = 'nonovl_'
        else:
            self.datasuf = ''
        if self.mixobs:
            self.datasuf = self.datasuf + 'on_mixed_impobs'
        else:
            self.datasuf = self.datasuf + 'on_imputed'
        # Model parameters:
        if model is None:
            self.model = 'observed_aux_18_on_mixed_impobs_QCUT'
        else:
            self.model = model
        self.tag = self.datasuf + '_MODEL_' + self.model
        self.chmmdir = os.environ['CHMMDIR'] + '/calls/' + self.model + '/'
        self.ecpref = self.chmmdir + 'ENH/' + self.model + \
            '_ENH_bin_' + self.datasuf + self.statesuf
        self.pcpref = self.chmmdir + 'PROM/' + self.model + \
            '_PROM_bin_' + self.datasuf + self.statesuf
        # TODO: Allow case where usemark is a list
        if self.usemark is not None:
            markdir = os.environ['IMPUTED_DIR'] + '/' + self.usemark + '/'
            self.hpref = markdir + self.usemark + '_all_bin_' + \
                self.datasuf + self.marksuf
            self.tag = self.tag + '_intersect_' + self.usemark
        # File and image prefixes
        self.filepref = os.environ['DBDIR'] + '/masterlist_matindices/' + \
            'matindices_' + self.tag
        self.figpref = os.environ['IMGDIR'] + '/masterlist_statistics/' + \
            'matindices_' + self.tag


    def main(self):
        self.load_data()
        self.calculate_margins()
        self.calculate_enhprom_ind()

    # Load in data and collapse matrix:
    def load_data(self):
        print("[STATUS] Loading data...")
        # Load data:
        self.eX = get_data(self.ecpref)
        self.pX = get_data(self.pcpref)
        # Load all names:
        self.enam = get_names(self.ecpref)
        self.pnam = get_names(self.pcpref)
        # Load histone matrix and intersect:
        if self.usemark is None:
            self.ehx = self.eX
            self.phx = self.pX
        else:
            self.hX = get_data(self.hpref)
            self.hnam = get_names(self.hpref)
            self.hnam = [re.sub(".*BSS","BSS", n).split("_")[0]
                         for n in self.hnam]
            # Reorder h matrix to match the other two:
            self.hord = np.array([np.where(np.array(self.hnam) == n)[0][0]
                                  for n in self.enam])
            self.hXord = self.hX[:, self.hord]
            # Turn hX into binary matrix, thresholded at 2:
            self.hXthresh = 1 * (self.hXord > 0)
            del(self.hX) # Remove the orig. histone matrix
            # Merge matrices:
            self.ehx = self.eX.multiply(self.hXthresh)
            self.phx = self.pX.multiply(self.hXthresh)
            print('[STATUS] Statistics on % kept dhs:')
            print('PROM x histone:', self.pX.nnz, "to", self.phx.nnz,
                  "is pct:", round(100 * self.phx.nnz / self.pX.nnz, 2))
            print('ENH x histone:', self.eX.nnz, "to", self.ehx.nnz,
                  "is pct:", round(100 * self.ehx.nnz / self.eX.nnz, 2))
            print('[STATUS] Statistics on % kept histone locations:')
            print('Histone x PROM:', self.hXthresh.nnz, "to", self.phx.nnz,
                  "is pct:", round(100 * self.phx.nnz / self.hXthresh.nnz, 2))
            print('Histone x ENH:', self.hXthresh.nnz, "to", self.ehx.nnz,
                  "is pct:", round(100 * self.ehx.nnz / self.hXthresh.nnz, 2))
            # Save matrices:
            if self.savematrices:
                save_X_names(self.ehx, self.enam, self.filepref +"_ENH_matrix",
                             savemtx=self.savemtx)
                save_X_names(self.phx, self.enam, self.filepref +"_PROM_matrix",
                             savemtx=self.savemtx)
                if self.saveraw:
                    # Also do this with strength of matrices:
                    self.ehx_raw = self.eX.multiply(self.hXord)
                    self.phx_raw = self.pX.multiply(self.hXord)
                    save_X_names(self.ehx_raw, self.enam,
                                 self.filepref + "raw_ENH_matrix",
                                 savemtx=self.savemtx)
                    save_X_names(self.phx_raw, self.enam,
                                 self.filepref + "raw_PROM_matrix",
                                 savemtx=self.savemtx)


    # Calculate margins, for evaluating enh/prom locations:
    def calculate_margins(self):
        print("[STATUS] Calculating margins")
        # Margins:
        self.le = np.array(np.sum(self.ehx, axis=0))[0]
        self.me = np.array(np.sum(self.ehx, axis=1)).T[0]
        self.lp = np.array(np.sum(self.phx, axis=0))[0]
        self.mp = np.array(np.sum(self.phx, axis=1)).T[0]
        # Locations where one or both are active
        self.df = pd.DataFrame({'ne': self.me, 'np': self.mp,
                               'sum': self.me + self.mp})
        # Total labeled with either enh/promoter:
        print(sum(self.df['sum'] > 0))
        # Locations where one is a
        self.indp = np.where(self.df['np'] > 0)[0]
        self.inde = np.where(self.df['ne'] > 0)[0]
        print(len(self.indp), "-", len(self.inde))
        if self.makeplots:
            self.plot_margins()


    def plot_margins(self):
        print("[STATUS] Plotting margins")
        # Plot one:
        fig = plt.figure(figsize=(14, 4))
        gs = gridspec.GridSpec(1, 2)
        # Enhancers subfigure:
        ax = plt.subplot(gs[0])
        u, c = np.unique(self.me, return_counts=True)
        plt.plot(np.cumsum(c[1:]))
        plt.ylabel('Cumulative number of enhancers')
        plt.xlabel('Number of epigenomes')
        # Promoters subfigure:
        ax = plt.subplot(gs[1])
        u, c = np.unique(self.mp, return_counts=True)
        plt.plot(np.cumsum(c[1:]))
        plt.ylabel('Cumulative number of promoters')
        plt.xlabel('Number of epigenomes')
        plt.tight_layout()
        fig = plt.gcf()
        fig.savefig(self.figpref + "_enhprom_margin_distr.png",
                    dpi=350, bbox_inches='tight')
        # Plot two:
        intp = [sum((self.df['ne'] <= i) & (self.df['np'] > 0))
                for i in range(21)]
        plt.bar(range(21), intp, color='darkgrey')
        plt.xlabel("At most X epigenomes labeled as enhancers")
        plt.ylabel("Number of regions")
        plt.grid(True, linestyle='dashed', color='darkgrey')
        fig = plt.gcf()
        fig.savefig(self.figpref + "_enhloc_margins.png",
                    dpi=350, bbox_inches='tight')


    # Choose enh and prom based on cutoffs of shared activity:
    def calculate_enhprom_ind(self):
        self.both = self.me + self.mp
        self.ind = np.where(self.both > 0)[0]
        self.fe = self.me[self.ind] / self.both[self.ind]
        self.fe = np.round(self.fe, 2)
        self.uf, self.cf = np.unique(self.fe, return_counts=True)
        self.ccf = np.cumsum(self.cf)
        # NOTE: We can calculate without H3K27ac as well by diff. run params
        self.enhind = self.ind[self.fe >= 0.75]
        self.promind = self.ind[self.fe <= 0.25]
        self.dyadicind= self.ind[(self.fe > 0.25) & (self.fe < 0.75)]
        print('Enhancers:', len(self.enhind))
        print('Promoters:', len(self.promind))
        print('Mixed:', len(self.ind) - (len(self.promind) + len(self.enhind)))
        print('Dyadic:', len(self.dyadicind))
        # Save enhancer, promoter, dyadic indices:
        print("[STATUS] Saving the enhancer and promoter indices")
        filesuf = '_masterlist_indices.cp.gz'
        AUX.save_pickle_gzip(self.filepref + '_ENH' + filesuf, self.enhind)
        AUX.save_pickle_gzip(self.filepref + '_PROM' + filesuf, self.promind)
        AUX.save_pickle_gzip(self.filepref + '_DYADIC' + filesuf, self.dyadicind)
        # Also save as tsv:
        tsvsuf = '_masterlist_indices.tsv'
        AUX.savelist(self.enhind, self.filepref + '_ENH' + tsvsuf)
        AUX.savelist(self.promind, self.filepref + '_PROM' + tsvsuf)
        AUX.savelist(self.dyadicind, self.filepref + '_DYADIC' + tsvsuf)
        if self.makeplots:
            self.plot_pct_each()
            self.plot_element_properties()
            self.plot_distr_enhprom()

    # Plots of which percentage enh/promoter:
    def plot_pct_each(self):
        print("[STATUS] Plotting pct each")
        # Get boxplot data:
        bi = self.both[self.ind]
        bpdata = []
        for x in tqdm(self.uf):
            idx = np.where(self.fe == x)[0]
            bpdata.append(bi[idx])
        # ------------------------------------------------------------
        f, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12,4))
        # Same data both:
        ax.bar(list(self.uf), list(self.cf), width=0.01)
        ax2.bar(list(self.uf), list(self.cf), width=0.01)
        # zoom-in / limit the view to different portions of the data
        t2 = np.sort(self.cf)[-2]
        ax.set_ylim(np.max(self.cf) * .9, np.max(self.cf) * 1.05)  # outliers only
        ax2.set_ylim(0, t2 * 1.2)  # most of the data
        # hide the spines between ax and ax2
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop='off')  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        # Format comma y-axis
        ax.get_yaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax2.get_yaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax.set_xlim(-0.01, 1.01)
        plt.xlabel('Fraction Annotated as Enhancer')
        add_diaglines(ax, ax2)
        f.text(0.02, 0.5, 'Number of Regions', va='center', rotation='vertical')
        ax.grid(True)
        ax2.grid(True)
        fig = plt.gcf()
        fig.savefig(self.figpref + "_frac_enhplot.png",
                    dpi=350, bbox_inches='tight')

        # ------------------------------------------------------------
        f, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12,6))
        # Boxplots of specificity:
        # pos = np.array(range(len(bpdata))) + 1
        pos = self.uf
        bp = ax.boxplot(bpdata, sym='k.', positions=pos, notch=1, widths=0.008)
        plt.setp(bp['whiskers'], color='k', linestyle='-')
        plt.setp(bp['fliers'], markersize=1.0)
        # this locator puts ticks at regular intervals
        # loc = mpl.ticker.MultipleLocator(base=0.2)
        # ax.xaxis.set_major_locator(loc)
        ax.xaxis.set_ticks(np.arange(0, 1.05, 0.2))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
        ax.set_ylabel('Number of Epigenomes')
        ax.set_xlim(-0.01, 1.01)
        ax.grid(True)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop='off')  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        #ax.set_xticks(np.arange(0, 1, 0.2))
        # Cumulative regions:
        ax2.bar(list(self.uf), list(self.ccf), width=0.01)
        ax2.set_ylim(0, self.ccf[-1] - .9 * self.cf[-1])  # most of the data
        ax2.get_yaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax2.grid(True)
        ax2.set_ylabel('Cumulative # Regions')
        ax2.set_xlabel('Fraction Annotated as Enhancer')
        fig = plt.gcf()
        fig.savefig(self.figpref + "_frac_enhplot2.png",
                    dpi=350, bbox_inches='tight')


    def plot_element_properties(self):
        print("[STATUS] Plotting element properties")
        # How close are promoters to each other?
        plt.figure(figsize=(6,4))
        dp = np.diff(self.promind)
        top = 20
        dp[dp > top] = top
        ud, cd = np.unique(dp, return_counts=True)
        ax = plt.gca()
        ax.bar(list(ud), list(cd), width=.8)
        ax.grid(True)
        ax.get_yaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x))))
        ax.set_ylabel('Number of Promoters')
        ax.set_xlabel('Distance (in masterlist rows) to next promoter')
        fig = plt.gcf()
        fig.savefig(self.figpref + "_prom_nearest_dist.png",
                    dpi=350, bbox_inches='tight')

        # How close are Enhancers to each other?
        plt.figure(figsize=(6,4))
        dp = np.diff(self.enhind)
        top = 20
        dp[dp > top] = top
        ud, cd = np.unique(dp, return_counts=True)
        ax = plt.gca()
        ax.bar(list(ud), list(cd), width=.8)
        ax.grid(True)
        ax.get_yaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x))))
        ax.set_ylabel('Number of Enhancers')
        ax.set_xlabel('Distance (in masterlist rows) to next enhancer')
        fig = plt.gcf()
        fig.savefig(self.figpref + "_enh_nearest_dist.png",
                    dpi=350, bbox_inches='tight')


    def plot_distr_enhprom(self):
        print("[STATUS] Recalculating margins")
        le = np.array(np.sum(self.ehx[self.enhind,:], axis=0))[0]
        lp = np.array(np.sum(self.phx[self.promind, :], axis=0))[0]
        # me = np.array(np.sum(self.ehx[self.enhind, :], axis=1)).T[0]
        # mp = np.array(np.sum(self.phx[self.promind, :], axis=1)).T[0]
        # Make plots:
        print("[STATUS] Plotting distribution: Number of enhancers, promoters")
        fig = plt.figure(figsize=(14, 4))
        gs = gridspec.GridSpec(1, 2)
        # Enhancers subfigure:
        ax = plt.subplot(gs[0])
        plt.hist(le, bins=40, color='darkgrey')
        plt.ylabel('Number of enhancers')
        plt.xlabel('Number of epigenomes')
        ax.get_yaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x))))
        # Promoters subfigure:
        ax = plt.subplot(gs[0])
        plt.hist(lp, bins=40, color='darkgrey')
        plt.ylabel('Number of enhancers')
        plt.xlabel('Number of epigenomes')
        ax.get_yaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x))))
        plt.tight_layout()
        fig = plt.gcf()
        fig.savefig(self.figpref + "_enhprom_final_margins.png",
                    dpi=350, bbox_inches='tight')


def get_data(pref):
    with gzip.open(pref + "_csr.cp.gz", 'rb') as f:
        X = pickle.load(f, encoding='latin1')
        return(X)


def get_names(pref):
    with gzip.open(pref + "_attr.cp.gz", 'rb') as f:
        attr = pickle.load(f, encoding='latin1')
    return(attr['names'])


def save_X_names(X, names, filepref, savemtx=False):
    data_file = filepref + "_csr.cp.gz"
    attr_file = filepref + "_attr.cp.gz"
    attr = {'names': names}
    print("[STATUS] Writing names to: " + attr_file)
    AUX.save_pickle_gzip(attr_file, attr)
    print("[STATUS] Writing X to: " + data_file)
    AUX.save_pickle_gzip(data_file, X)
    # With option of saving in more general format:
    if savemtx:
        data_file = filepref + ".mtx"
        attr_file = filepref + "_names.tsv"
        print("[STATUS] Also saving to mtx file: " + data_file)
        AUX.savelist(names, attr_file)
        mmwrite(data_file, X)


# For broken plot:
def add_diaglines(ax, ax2, d=0.015):
    #d = .015  # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal


if __name__ == "__main__":
    fire.Fire(calculate_functional_idx)
