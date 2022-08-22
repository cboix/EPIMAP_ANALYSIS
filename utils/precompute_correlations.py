#!/usr/bin/env python3
# -----------------------------------------------------------------------
# Pre-compute the correlations between elements and genes:
# MODE 1. Precompute correlation for any element within 1,010kb of a gene
# MODE 2. Precompute N correlations with random genes for each enhancer
# -----------------------------------------------------------------------
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
from scipy.sparse import csr_matrix

import socket

domain = socket.getfqdn()
import matplotlib as mpl
if 'broadinstitute.org' in domain:
    mpl.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns

from sklearn.metrics import average_precision_score, precision_recall_curve

import parse_config
import auxfunc_masterlist as AUX


class precompute_corr(object):
    def __init__(self, exprfile, markfile, tssfile, mlocfile, outprefix,
                 indfile=None, pairsfile=None, window=1.010**6,
                 numrand=25, random=False):
        self.exprfile = exprfile
        self.markfile = markfile
        self.tssfile = tssfile
        self.mlocfile = mlocfile
        self.outprefix = outprefix
        self.indfile = indfile  # For keep specific ind
        self.pairsfile = pairsfile # For possible pairs
        self.window = window
        self.datadir = os.environ['DBDIR']
        self.random = random
        self.numrand = numrand

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
            all_genes = self.tssdf.gene.to_numpy()
            NGENES = len(all_genes)
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
                    if self.random:
                        # Sample self.numrand genes per enhancer and append:
                        for mi in range(NM):
                            mtmp = markloc[mi]
                            erand = np.random.choice(NGENES, self.numrand)
                            for ei in erand:
                                gtmp = all_genes[ei]
                                pdlist.append([gtmp, mtmp])
                    else:
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
                    if self.random:
                        dflist[chrom] = pd.DataFrame(pdlist, columns=['gene','mid'])
                        dflist[chrom]['chrom'] = chrom
                    else:
                        dflist[chrom] = pd.DataFrame(pdlist, columns=['tss','mid'])
                        dflist[chrom]['chrom'] = chrom
            df_aslist = [v for k,v in dflist.items()]
            self.pairsdf = pd.concat(df_aslist)
            print(self.pairsdf.head())
            self.mlocdf['mind'] = self.mlocdf.index
            print(self.mlocdf.head())
            self.pairsdf = self.pairsdf.merge(self.mlocdf)
            print(self.pairsdf.head())
            # NOTE: Don't merge with genes if random, arent on same chr.
            if not self.random:
                self.pairsdf = self.pairsdf.merge(self.tssdf)
                print(self.pairsdf.head())
            print(self.pairsdf.shape)
            if self.pairsfile is not None:
                print("[STATUS] Saving to file:", self.pairsfile)
                AUX.save_pickle_gzip(self.pairsfile, self.pairsdf)

    # Calculate the correlation for each pair
    # TODO: Do so for all marks and save as a feature matrix
    def calc_correlation(self):
        if self.random:
            corrfile = self.outprefix + "_random" + "_corr.tsv.gz"
        else:
            corrfile = self.outprefix + "_corr.tsv.gz"
        if os.path.exists(corrfile):
            print('[STATUS] Loading corr from file:', corrfile)
            self.allcorr = pd.read_csv(corrfile, header=None,
                                       sep="\t").to_numpy()
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
                # Conv if class sparse:
                if type(subX) == csr_matrix:
                    subX = subX.toarray()
                sgenes = sub_pdf.gene.to_numpy()
                subE = self.exprwide.loc[sgenes]
                subE = np.array(subE)
                corr_vec = np.array(pearson_mat(subX, subE))
                self.allcorr[indlist] = corr_vec
            np.savetxt(corrfile, self.allcorr)

    def main(self):
        self.load_data()
        self.setup_pairsdf()  # Sets up all testable links (reusable)
        self.calc_correlation()


def fmtcm(x):
    return("{:,}".format(int(x)))


# TODO: Also allow pre-compute KL.
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


if __name__ == "__main__":
    fire.Fire(precompute_corr)
