#!/usr/bin/python
# ------------------------------------------------------------------
# Utility to return genome-wide statistics imputed recover of peaks:
# ------------------------------------------------------------------
import re
import os
import sys
import fire
import gzip
import numpy as np
import pandas as pd
import sklearn.metrics

import parse_config

# TODO: Match imputed and observed data,


class calculate_gw_recovery(object):
    def __init__(self, mark, infofile, dir, outprefix, chrom=None):
        # Arguments:
        self.mark = mark
        self.infofile = infofile  # Info table - what to run
        self.dir = dir  # Directory where files located
        self.outprefix = outprefix  # Output prefix
        self.chrom = chrom
        if self.chrom is None:
            self.chrlist = ['chr' + str(i+1) for i in range(22)] + ['chrX']
        else:
            if type(chrom) == str:
                self.chrlist = [chrom]
            else:
                self.chrlist = chrom
        # Files and directories:
        self.datadir = os.environ['DBDIR']

    def main(self):
        self.process_infofile()
        for chrom in self.chrlist:
            self.get_filelist(chrom)
            self.init_matrix()
            self.read_data()
            self.calculate_recovery()
            self.calculate_recovery_bestobs()
            self.write_statistics(chrom)

    def process_infofile(self):
        self.info= pd.read_csv(
            self.infofile, skiprows=0, header=None, sep="\t",
            names=['id','mark','prefix'])
        self.info = self.info[self.info.mark == self.mark]
        self.info['type'] = ['Observed' if re.match('FINAL', v)
                             else 'Imputed' for v in self.info.prefix]
        print(self.info.head())
        print(self.info.shape)
        # --------------------------------
        # Match the ones with the same id:
        u, c = np.unique(self.info['id'], return_counts=True)
        self.matched_names = u[c == 2]
        self.nmatched = len(self.matched_names)
        print('nmatched', self.nmatched)
        print(self.matched_names[0:5])
        # -----------------------------------
        # Reduce to only observed or matched:
        keepind = ((self.info.id.isin(self.matched_names)) |
                   (self.info['type'] == 'Observed'))
        self.info = self.info[keepind]
        self.info = self.info.reset_index(drop=True)
        print('Reduced info to:', self.info.shape)
        # ----------------------------------
        # Get indices of observed data only:
        self.obs_inds = self.info.index[self.info['type'] == 'Observed'].tolist()
        print('nobs', len(self.obs_inds))
        print((self.info['id'] == self.matched_names[1]) & (self.info['type'] == 'Observed'))
        # -------------------------------------------------
        # Get the observed and imputed indices for matched:
        self.mobs = [np.where(((self.info['id'] == self.matched_names[i]) &
            (self.info['type'] == 'Observed')).tolist())[0][0]
                     for i in range(self.nmatched)]
        self.mimp = [np.where(((self.info['id'] == self.matched_names[i]) &
            (self.info['type'] == 'Imputed')).tolist())[0][0]
                     for i in range(self.nmatched)]
        print(self.mobs)
        print(self.mimp)
        print(self.info.loc[self.mobs[0:5]])
        print(self.info.loc[self.mimp[0:5]])

    def get_filelist(self, chrom):
        self.filelist = [chrom + "_" + p + ".wig.gz" for p in self.info.prefix]
        self.nfile = len(self.filelist)
        print("[STATUS] Will read in " + str(self.nfile) + " files.")

    def init_matrix(self):
        fname = self.filelist[1]
        x = pd.read_csv(self.dir + "/" + fname, skiprows=2, header=None).to_numpy()
        self.nrow = len(x)
        self.X = np.zeros((self.nrow, self.nfile))

    def read_data(self):
        for i, fname in enumerate(self.filelist):
            print(str(i) + " " + fname)
            x = pd.read_csv(self.dir + "/" + fname, skiprows=2, header=None).to_numpy()
            self.X[:,i] = x.T[0]
        print(self.X[0:10,])
        print(self.X.shape)
        # Calculate the mean of observed:
        self.mX = np.mean(self.X[:, np.array(self.mobs)], axis=1)

    def calculate_recovery(self):
        self.meanZ = []
        self.ccf = []
        self.auroc2 = []
        self.aurocq = []
        self.ap2 = []
        self.apq = []
        for i in range(self.nmatched):
            oind = self.mobs[i]
            iind = self.mimp[i]
            samp = self.matched_names[i]
            print(i, samp, oind, iind)
            # Corr, auroc, auprc, etc.
            oX = self.X[:,oind]
            iX = self.X[:,iind]
            oZ = 1 * (oX > 2)
            # Recovery of the top 1% of peaks:
            thresh = np.quantile(oX, 0.99)
            oZq = 1 * (oX > thresh)
            self.meanZ.append(np.mean(oZ))
            print(np.mean(oZq))
            self.ccf.append(np.corrcoef(oX, iX)[0,1])
            # AUCs:
            fpr, tpr, thresholds = sklearn.metrics.roc_curve(
                y_true=oZ, y_score=iX, pos_label=1)
            self.auroc2.append(sklearn.metrics.auc(fpr, tpr))
            fpr, tpr, thresholds = sklearn.metrics.roc_curve(
                y_true=oZq, y_score=iX, pos_label=1)
            self.aurocq.append(sklearn.metrics.auc(fpr, tpr))
            self.ap2.append(sklearn.metrics.average_precision_score(oZ, iX))
            self.apq.append(sklearn.metrics.average_precision_score(oZq, iX))

    # Calculate recovery metrics for the BEST matched neighbor:
    def calculate_recovery_bestobs(self):
        self.mean_auroc2 = []
        self.mean_aurocq = []
        self.mean_ap2 = []
        self.mean_apq = []
        self.near_auroc2 = []
        self.near_aurocq = []
        self.near_ap2 = []
        self.near_apq = []
        self.near_samp = []
        self.corrobs = np.corrcoef(self.X[:, self.obs_inds].T)
        print("corrobs.shape", self.corrobs.shape)
        for i in range(self.nmatched):
            oind = self.mobs[i]
            samp = self.matched_names[i]
            print(i, samp, oind)
            oX = self.X[:,oind]
            oZ = 1 * (oX > 2)
            thresh = np.quantile(oX, 0.99)
            oZq = 1 * (oX > thresh)
            # Mean metrics:
            fpr, tpr, thresholds = sklearn.metrics.roc_curve(
                y_true=oZ, y_score=self.mX, pos_label=1)
            self.mean_auroc2.append(sklearn.metrics.auc(fpr, tpr))
            fpr, tpr, thresholds = sklearn.metrics.roc_curve(
                y_true=oZq, y_score=self.mX, pos_label=1)
            self.mean_aurocq.append(sklearn.metrics.auc(fpr, tpr))
            self.mean_ap2.append(sklearn.metrics.average_precision_score(
                oZ, self.mX))
            self.mean_apq.append(sklearn.metrics.average_precision_score(
                oZq, self.mX))
            # Reduce the tested corr top ten by corr:
            oi = np.array(np.where(self.obs_inds == oind)[0])
            print(i, oind, oi)
            corrvec = np.array(self.corrobs[oi,:].tolist()[0])
            thresh = np.sort(corrvec)[-10]
            tind = np.where(corrvec >= thresh)[0]
            testids = [self.obs_inds[t] for t in tind]
            print('closest:', testids)
            # Nearest metrics:
            top_apq = 0
            top_samp = ''
            for j in testids:
                jsamp = self.info['id'].tolist()[j]
                if jsamp != samp:
                    jX = self.X[:, j]
                    apj = sklearn.metrics.average_precision_score(oZq, jX)
                    if apj > top_apq:
                        print(j, jsamp, apj)
                        top_apq = apj
                        top_samp = jsamp
                        top_ap2 = sklearn.metrics.average_precision_score(
                            oZ, jX)
                        fpr, tpr, thresholds = sklearn.metrics.roc_curve(
                            y_true=oZq, y_score=jX, pos_label=1)
                        top_aurocq = sklearn.metrics.auc(fpr, tpr)
                        fpr, tpr, thresholds = sklearn.metrics.roc_curve(
                            y_true=oZ, y_score=jX, pos_label=1)
                        top_auroc2 = sklearn.metrics.auc(fpr, tpr)
            # Add these nearest samples
            self.near_samp.append(top_samp)
            self.near_apq.append(top_apq)
            self.near_ap2.append(top_ap2)
            self.near_aurocq.append(top_aurocq)
            self.near_auroc2.append(top_auroc2)
            # Print stats as we go:
            print('Mean:', samp, 'mean', self.mean_aurocq[-1], self.mean_ap2[-1])
            print('Top:', samp, top_samp, top_aurocq, top_apq)

    def write_statistics(self, chrom):
        # Make dictionary:
        self.statsdf = pd.DataFrame({
            'id': self.matched_names, 'meanZ': self.meanZ, 'cor': self.ccf,
            'auroc2': self.auroc2, 'aurocq': self.aurocq,
            'ap2': self.ap2, 'apq': self.apq,
            'mean_auroc2': self.mean_auroc2, 'mean_aurocq': self.mean_aurocq,
            'mean_ap2': self.mean_ap2, 'mean_apq': self.mean_apq,
            'near_id': self.near_samp,
            'near_auroc2': self.near_auroc2, 'near_aurocq': self.near_aurocq,
            'near_ap2': self.near_ap2, 'near_apq': self.near_apq})
        self.statsdf.to_csv(self.outprefix + "_" + chrom + "_stats.tsv.gz",
                            sep="\t", index=False)


def pearson_mat(A, B):
    # Row-wise mean of input arrays & subtract from input arrays themselves
    A_mA = A - A.mean(1)[:,None]
    B_mB = B - B.mean(1)[:,None]
    # Sum of squares across rows
    ssA = (A_mA**2).sum(1);
    ssB = (B_mB**2).sum(1);
    out = (A_mA * B_mB).sum(1) / np.sqrt(ssA * ssB)
    return(out)

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


if __name__ == "__main__":
    fire.Fire(calculate_gw_recovery)

