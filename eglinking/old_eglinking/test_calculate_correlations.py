#!/usr/bin/env python3
# Calculate + plot correlations between local elements:
# TODO: COMPARE TO THE ROADMAP LINKS
import os
import re
import sys
import time
import gzip
import fire

import itertools
import pickle
import numpy as np
import pandas as pd

import socket

domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')
# from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import seaborn as sns

from sklearn.metrics import average_precision_score, precision_recall_curve

# import json
# from sklearn.linear_model import LogisticRegression
# from scipy.stats import pearsonr

import parse_config
import auxfunc_masterlist as AUX



REGULARIZATION=1.0
PENALTY="l2"


class calc_corr(object):
    def __init__(self, exprfile, markfile, tssfile, mlocfile, outprefix,
                 indfile=None, pairsfile=None, window=10**6):
        self.exprfile = exprfile
        self.markfile = markfile
        self.tssfile = tssfile
        self.mlocfile = mlocfile
        self.outprefix = outprefix
        self.indfile = indfile  # For keep specific ind
        self.pairsfile = pairsfile # For possible pairs
        self.window = window
        self.datadir = os.environ['DBDIR']

    def read_tss(self):
        print("[STATUS] Loading the TSS data:", self.tssfile)
        self.tssdf = pd.read_csv(self.tssfile, header=None, sep='\t',
                                 names=['gene', 'chrom','tss','strand'])
        self.chromlist = pd.unique(self.tssdf['chrom'])
        print(self.chromlist)

    # Read in the expression data:
    def read_expr(self):
        self.exprpref = re.sub(".mtx.gz","", self.exprfile)
        self.exprwidefile = self.exprpref + ".wide.cp.gz"
        if not os.path.exists(self.exprwidefile):
            print("[STATUS] Loading the expression data from MTX")
            self.expr_id = []
            self.expr_gene = []
            self.expr_val = []
            self.exprdf = pd.read_csv(self.exprfile, sep='\t')
            print(self.exprdf.shape)
            print(self.exprdf.head())
            # Pivot into wide table:
            self.exprwide = self.exprdf.pivot(index='gene', columns='id', values='log2fpkm')
            print(self.exprwide.head())
            print(self.exprwide.shape)
            # Save matrix:
            print("[STATUS] Saving as pickled objects")
            AUX.save_pickle_gzip(self.exprwidefile, self.exprwide)
        else:
            print("[STATUS] Loading the expression data from pickle files")
            self.exprwide = AUX.load_pickle_gzip(self.exprwidefile)

    def read_markloc(self):
        # Read in file:
        self.mlocdf = pd.read_csv(self.mlocfile, sep="\t", header=None,
                                  names=['chrom','start','end','name'])
        self.mlocdf['mid'] = (self.mlocdf['start'] + self.mlocdf['end']) / 2

    def read_mark(self):
        self.markpref = re.sub("_csr.cp.gz","", self.markfile)
        print("[STATUS] Loading the mark file:", self.markfile)
        self.X = AUX.load_pickle_gzip(self.markfile)
        print("[STATUS] Loaded file of size:", self.X.shape)
        attr = AUX.load_pickle_gzip(self.markpref + "_attr.cp.gz")
        self.xnames = [re.sub("_.*\n","", re.sub(".*BSS","BSS", nam)) for nam in attr['names']]

    def load_data(self):
        # Read expr data:
        self.read_expr()
        self.read_tss()
        # Read mark data:
        self.read_mark()
        self.read_markloc()
        # Match names:
        self.exprnames = list(self.exprwide.columns)
        print(self.exprnames[0:5])
        print(self.xnames[0:5])
        self.kept_names = [nam for nam in self.xnames if nam in self.exprnames]
        print("[STATUS] Keeping", len(self.kept_names), "samples:",
              self.kept_names[0:20], "etc...")
        self.exprwide = self.exprwide[self.kept_names]
        print(self.exprwide.head())
        self.kept_idx = [i for i, nam in enumerate(self.xnames) if nam in self.kept_names]
        print(self.kept_idx[0:5])
        self.X = self.X[:, self.kept_idx]
        print("[STATUS] New X shape:", self.X.shape)
        # Slice data by indices:
        if self.indfile is not None:
            print("[STATUS] Reducing mark file to indices in:", self.indfile)
            self.markind = AUX.load_pickle_gzip(self.indfile)
            self.X = self.X[self.markind,:]
            print(self.X.shape)
            self.mlocdf = self.mlocdf.iloc[self.markind]
            self.mlocdf.reset_index(inplace=True, drop=True)
            print(self.mlocdf.shape)

    def setup_pairsdf(self):
        if self.pairsfile is not None and os.path.exists(self.pairsfile):
            print("[STATUS] Loading pairsdf from file:", self.pairsfile)
            self.pairsdf = AUX.load_pickle_gzip(self.pairsfile)
        else:
            # Find all possible within self.window...
            dflist = {}
            for chrom in self.chromlist:
                print(chrom)
                eind = np.where(self.tssdf.chrom == chrom)[0]
                exprloc = self.tssdf.tss[eind].to_numpy()
                #print(eind[0:5], exprloc[0:5])
                mind = np.where(self.mlocdf.chrom == chrom)[0]
                #print(self.mlocdf.head())
                markloc = self.mlocdf.mid[mind].to_numpy()
                use_iter = False
                if use_iter:
                    #print(mind[0:5], markloc[0:5])
                    params = itertools.product(exprloc, markloc)
                    comp = (lambda em: np.abs(em[0] - em[1]) < self.window)
                    filtered_params = itertools.filterfalse(comp, params)
                    tmpmat  = list(filtered_params)
                    dflist[chrom] = pd.DataFrame(tmpmat, columns=['eloc','mloc'])
                    dflist[chrom]['chrom'] = chrom
                else:
                    NM = len(markloc)
                    NE = len(exprloc)
                    # Ensure sorted:
                    markloc = np.sort(markloc)
                    exprloc = np.sort(exprloc)
                    # For each gene in the chromosome:
                    firstenh = 0
                    pdlist = []
                    for ei in range(NE):
                        foundany = 0
                        etmp = exprloc[ei]
                        for mi in range(firstenh, NM):
                            mtmp = markloc[mi]
                            dist = mtmp - etmp
                            if np.abs(dist) <= self.window:
                                pdlist.append([etmp, mtmp])
                                # Flag first enh to start here next time
                                if not foundany:
                                    foundany = 1
                                    firstenh = mi
                            # If past window, break:
                            elif dist > self.window:
                                break
                    dflist[chrom] = pd.DataFrame(pdlist, columns=['tss','mid'])
                    dflist[chrom]['chrom'] = chrom
            df_aslist = [v for k,v in dflist.items()]
            self.pairsdf = pd.concat(df_aslist)
            print(self.pairsdf.head())
            self.mlocdf['mind'] = self.mlocdf.index
            print(self.mlocdf.head())
            self.pairsdf = self.pairsdf.merge(self.mlocdf)
            print(self.pairsdf.head())
            self.pairsdf = self.pairsdf.merge(self.tssdf)
            print(self.pairsdf.head())
            if self.pairsfile is not None:
                print("[STATUS] Saving to file:", self.pairsfile)
                AUX.save_pickle_gzip(self.pairsfile, self.pairsdf)


    # Calculate the correlation for each pair
    # TODO: Do so for all marks and save as a feature matrix
    def calc_correlation(self):
        corrfile = self.outprefix + "_corr.tsv"
        if os.path.exists(corrfile):
            print('[STATUS] Loading corr from file:', corrfile)
            self.allcorr = pd.read_csv(corrfile, header=None, sep="\t").to_numpy()
        else:
            npairs = self.pairsdf.shape[0]
            print('[STATUS] Calc. corr for', npairs, 'pairs')
            chunksize = int(1e5)
            nchunk = int(np.floor(npairs / chunksize)) + 1
            self.allcorr = np.zeros(npairs)
            print(nchunk)
            for i in range(nchunk):
                print(i)
                indlist = list(range(i * chunksize,
                                     np.min([(i+1) * chunksize, npairs])))
                sub_pdf = self.pairsdf.iloc[indlist]
                sloc = sub_pdf.mind.to_numpy()
                subX = self.X[sloc,:]
                subX = subX.toarray()
                sgenes = sub_pdf.gene.to_numpy()
                subE = self.exprwide.loc[sgenes]
                subE = np.array(subE)
                corr_vec = np.array(pearson_mat(subX, subE))
                self.allcorr[indlist] = corr_vec
            np.savetxt(corrfile, self.allcorr)

    def main(self):
        self.load_data()
        self.setup_pairsdf()  # Sets up all testable links:
        self.calc_correlation()
        self.load_pchic()

    # TODO: Calc corr with closest:


    # Load PromoterCapture HiC links, match to pairsdf
    def load_pchic(self):
        # TODO REMOVE
        p
        self.lddir = self.datadir + '/linking_data/'
        self.pcfile = self.lddir + \
            'ActivePromoterEnhancerLinks_JavierrePChIC.tsv.gz'
        self.pcdf = pd.read_csv(self.pcfile, sep="\t")
        self.pcdf = self.pcdf.loc[:, ['baitID', 'oeID']].drop_duplicates()
        # Read the overlaps of
        self.btdf = pd.read_csv(self.lddir + 'pchic_bait_gene.tsv', sep="\t")
        # TODO: OEINT is for the OVERLAPPING DATASET NOT NON-OVERLAPPING
        self.oedf = pd.read_csv(self.lddir + 'oeint.tsv', sep="\t")
        # Merge back to matches:
        self.pcdf = self.pcdf.merge(self.btdf)
        self.pcdf = self.pcdf.merge(self.oedf)
        # Keep only the attr that will merge into pairsdf
        self.pcdf = self.pcdf.loc[:, ['name', 'gene']]
        # self.pcdf = self.pcdf.drop_duplicates()
        print("[STATUS] Unique promoter - enhancer PCHiC interactions: " +
              fmtcm(len(self.pcdf)))
        # Merge with the test df:
        self.pcdf['ispc'] = 1
        self.pairsdf['ord'] = self.pairsdf.index
        print(self.pairsdf.shape)
        self.allcorr[np.isnan(self.allcorr)] = 0
        self.pairsdf['pcorr'] = self.allcorr
        self.pairsdf = self.pairsdf.merge(self.pcdf, how='left')
        self.pairsdf.ispc.fillna(0, inplace=True)
        print(self.pairsdf.shape)
        # Plot the links distributions:
        self.pairsdf['dist'] = self.pairsdf.mid - self.pairsdf.tss
        self.plot_dist_histogram(self.pairsdf.dist.to_numpy(),
                                 filename='dist_all.png',
                                 title='All tested links')
        pcdist = self.pairsdf.dist[self.pairsdf.ispc == 1]
        self.plot_dist_histogram(pcdist.to_numpy(),
                                 filename='dist_pchic.png',
                                 title='Links matching PC-HiC')

        # Plot the distibution of correlation values:
        self.plot_dist_histogram(self.pairsdf.pcorr.to_numpy(), isdist=False,
                                 filename='corr_all.png',
                                 title='Links matching PC-HiC')
        pc_corr = self.pairsdf.pcorr[self.pairsdf.ispc == 1]
        self.plot_dist_histogram(pc_corr.to_numpy(), isdist=False,
                                 filename='corr_pchic.png',
                                 title='Links matching PC-HiC')

        # TODO: Reorder if not ordered.
        print(self.pairsdf.head())
        print("[STATUS] Merge with test links: " +
              fmtcm(np.sum(self.pairsdf.ispc)) + " PCHiC links in " +
              fmtcm(len(self.pairsdf)) + " total possible links")
        y_true = self.pairsdf.ispc.to_numpy().astype(int)
        print(y_true.shape)
        pcorr = self.pairsdf.pcorr.to_numpy()
        # Score with correlation:
        score1_corr = (pcorr - np.min(pcorr)).flatten()
        score1_corr = score1_corr / np.max(score1_corr)
        print(score1_corr.shape)
        ap1 = average_precision_score(y_true, score1_corr)
        print(ap1)
        self.plot_auprc(y_true, score1_corr, filename='auprc_test1.png', title='Corr Only')
        # Score with distance
        score2_dist = np.exp(- np.abs(self.pairsdf.mid - self.pairsdf.tss) / 20000).to_numpy().flatten()
        score2_dist = score2_dist / np.max(score2_dist)
        print(score2_dist.shape)
        ap2 = average_precision_score(y_true, score2_dist)
        print(ap2)
        self.plot_auprc(y_true, score2_dist, filename='auprc_test2.png', title='Dist Only')
        # Score with both:
        score3_both = score1_corr * score2_dist
        score3_both = score3_both / np.max(score3_both)
        print(score3_both.shape)
        # Calculate precision recall curve and AP scores:
        ap3 = average_precision_score(y_true, score3_both)
        print("[STATUS] AP:", [ap1, ap2, ap3])
        # precision, recall, thresholds = precision_recall_curve(y_true, score3_both)
        # PLOT:
        self.plot_auprc(y_true, score3_both, filename='auprc_test3.png', title='Both Corr and Dist')

    def plot_auprc(self, y_true, y_test, filename, title=None):
        precision, recall, thresholds = precision_recall_curve(y_true, y_test)
        apval = average_precision_score(y_true, y_test)
        sns.set(font_scale=1.1)
        # Plot current trace:
        fig = plt.figure(figsize=(9, 6))
        ax = plt.gca()
        # Plot % change:
        if title is not None:
            ax.set_title("Precision-Recall Curve for " + str(title) + " (" + str(apval) + ")")
        ax.set_facecolor('white')
        plt.plot(recall, precision, lw=2)
        plt.ylabel('Precision')
        plt.xlabel('Precision')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        plt.tight_layout()
        fig = plt.gcf()
        fig.savefig(filename, dpi=350, bbox_inches='tight')

    def plot_dist_histogram(self, dist, filename, title=None, isdist=True):
        fig = plt.figure(figsize=(9, 5))
        ax = plt.gca()
        # Enhancers subfigure:
        plt.hist(dist, bins=40, color='darkgrey')
        plt.ylabel('Number of links')
        ax.get_yaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        # Label as either dist or corr:
        if isdist:
            plt.xlabel('Distance between link endpoints')
            ax.get_xaxis().set_major_formatter(
                mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        else:
            plt.xlabel('Correlation between mark and expression in link')
        plt.tight_layout()
        fig = plt.gcf()
        fig.savefig(filename, dpi=350, bbox_inches='tight')


def fmtcm(x):
    return("{:,}".format(int(x)))


def kullback_leibler(A, B):
    total = 0
    for a, b in zip(A, B):
        if b != 0 and a != 0: # if a=0 then total is unchanged
            lg = np.log(a/b)
            total += a * lg
    return total


def pearson_mat(A, B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:,None]
    B_mB = B - B.mean(1)[:,None]
    # Sum of squares across rows
    ssA = (A_mA**2).sum(1);
    ssB = (B_mB**2).sum(1);
    out = (A_mA * B_mB).sum(1) / np.sqrt(ssA * ssB)
    return(out)


# pos_corr_vec = np.array(pearson_mat(xfull, expr_list))
# neg_corr_vec = np.array(pearson_mat(xfull, rand_list))
# pos_corr_vec[np.isnan(pos_corr_vec)] = 0
# neg_corr_vec[np.isnan(neg_corr_vec)] = 0
# pos_data[:,j] = pos_corr_vec
# neg_data[:,j] = neg_corr_vec


if __name__ == "__main__":
    fire.Fire(calc_corr)
