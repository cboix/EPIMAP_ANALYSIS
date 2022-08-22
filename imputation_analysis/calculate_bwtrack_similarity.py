#!/usr/bin/python
# ------------------------------------------------------------------------
# Utility alternative to gwcorr from ChromImpute
# ------------------------------------------------------------------------
import re
import os
import sys
import fire
import gzip
import numpy as np
import pandas as pd

import parse_config


# TODO: CALC bw sim 1 track to all.
# Currently no changes - just copied over calc mark var.
class calculate_bw_similarity(object):
    def __init__(self, task, infofile, dir, outprefix,
                 dataset=None, chrom=None):
        # Arguments:
        self.task = task
        self.infofile = infofile  # Info table - what to run
        self.dir = dir  # Directory where files located
        self.outprefix = outprefix  # Output prefix
        self.chrom = chrom
        self.dataset = dataset
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
        self.info= pd.read_csv(
            self.infofile, skiprows=0, header=None, sep="\t",
            names=['id','mark','prefix'])
        self.info = self.info[self.info.mark == self.mark]
        if self.dataset == 'observed':
            ind = [i for i, v in enumerate(self.info.prefix) if re.match('FINAL', v)]
            self.info = self.info.iloc[np.array(ind)]
        elif self.dataset == 'imputed':
            ind = [i for i, v in enumerate(self.info.prefix) if re.match('impute', v)]
            self.info = self.info.iloc[np.array(ind)]
        for chrom in self.chrlist:
            self.get_filelist(chrom)
            self.init_matrix()
            self.read_data()
            self.calculate_statistics()
            self.write_statistics(chrom)

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

    def calculate_statistics(self):
        self.mX = np.mean(self.X, axis=1)
        self.sX = np.std(self.X, axis=1)
        self.cX = self.sX / self.mX
        self.qX = np.sum(self.X >= 2, axis=1)  # Num peaks

    def write_statistics(self, chrom):
        gzipped_write(self.mX, self.outprefix + "_" + chrom + "_mean.tsv.gz")
        gzipped_write(self.sX, self.outprefix + "_" + chrom + "_std.tsv.gz")
        gzipped_write(self.cX, self.outprefix + "_" + chrom + "_coeffv.tsv.gz")
        gzipped_write(self.qX, self.outprefix + "_" + chrom + "_npeak.tsv.gz")

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
    fire.Fire(calculate_bw_similarity)

