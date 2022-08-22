#!/usr/bin/pyton
# ------------------------------------------------------------------------
# Script to return the average track for each group + datatype + mark set:
# ------------------------------------------------------------------------
# Use: all_released_tracks_withgroup.tsv
import re
import os
import sys
import fire
import gzip
import numpy as np
import pandas as pd

import parse_config


class calculate_average_bw(object):
    def __init__(self, mark, infofile, dir, outprefix,
                 dataset='observed', group=None, chrom=None):
        # Arguments:
        self.mark = mark
        self.infofile = infofile  # Info table - what to run
        self.dir = dir  # Directory where files located
        self.outprefix = outprefix  # Output prefix
        self.chrom = chrom
        self.dataset = dataset
        self.group = group
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
            names=['id','mark','prefix', 'group'])
        self.info = self.info[self.info.mark == self.mark]
        # Pick only observed/imputed, etc:
        if self.dataset == 'observed':
            ind = [i for i, v in enumerate(self.info.prefix)
                   if re.match('FINAL', v)]
            self.info = self.info.iloc[np.array(ind)]
        elif self.dataset == 'imputed':
            ind = [i for i, v in enumerate(self.info.prefix)
                   if re.match('impute', v)]
            self.info = self.info.iloc[np.array(ind)]
        # Subset by group:
        if self.group is not None:
            self.info = self.info.loc[self.info.group == self.group]
        print(self.info.shape)
        print(self.info.head())
        for chrom in self.chrlist:
            self.get_filelist(chrom)
            self.init_matrix()
            self.read_data()
            # TODO: Here, add calculate average tracks for each group
            # separately.
            self.calculate_averages()
            self.write_averages(chrom)

    def get_filelist(self, chrom):
        self.filelist = [chrom + "_" + p + ".wig.gz" for p in self.info.prefix]
        self.nfile = len(self.filelist)
        print("[STATUS] Will read in " + str(self.nfile) + " files.")

    def init_matrix(self):
        fname = self.filelist[0]
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

    def calculate_averages(self):
        if self.group is not None:
            self.mX = np.mean(self.X, axis=1)
            print(self.mX.shape)
        else:
            self.mdict = {}
            # Calculate a per group average:
            for group in pd.unique(self.info.group):
                print(group)
                ind = [i for i, v in enumerate(self.info.group) if v == group]
                print(ind)
                self.mdict[group] = np.mean(self.X[:,ind], axis=1)
                print(self.mdict[group].shape)

    def write_averages(self, chrom):
        if self.group is not None:
            groupstr = re.sub(" & ", "and", self.group)
            gzipped_write(self.mX, self.outprefix + "_" +
                          groupstr + "_" + chrom + "_mean.tsv.gz")
        else:
            # Calculate a per group average:
            for group in pd.unique(self.info.group):
                print(group)
                groupstr = re.sub(" & ", "and", group)
                gzipped_write(self.mdict[group], self.outprefix + "_" +
                              groupstr + "_" + chrom + "_mean.tsv.gz")


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
    fire.Fire(calculate_average_bw)

