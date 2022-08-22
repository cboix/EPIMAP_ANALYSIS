#!/usr/bin/python
# --------------------------------------------------
# Count number of overlaps for each motif + shuffles
# In a set of regions and a background set
# NOTE: STILL WAYYYY TOO SLOW.
# --------------------------------------------------
import os
import re
import gzip
import fire
import time
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

rng = np.random


class count_motif_instances(object):
    def __init__(self, infile, bkgfile, outfile, motif_list, motif_pref):
        # Arguments:
        self.infile = infile
        self.bkgfile = bkgfile
        self.outfile = outfile
        # Motif parameters
        self.motif_list = motif_list
        self.motif_pref = motif_pref
        self.chrlist = ['chr' + str(i+1) for i in range(22)] + ['chrX']

    def main(self):
        self.read_data()
        self.read_motif_info()
        self.count_motifs()
        self.write_counts()

    def read_data(self):
        t1 = time.time()
        print("[STATUS] Reading input file")
        self.indf = pd.read_csv(self.infile, skiprows=0, header=None,
                                sep="\t", names=['chr','start','end','n'])
        print(self.indf.shape)
        print("[STATUS] Reading background file")
        self.bkgdf = pd.read_csv(self.bkgfile, skiprows=0, header=None,
                                sep="\t", names=['chr','start','end','n'])
        print(self.bkgdf.shape)
        self.ist = {}
        self.iend = {}
        self.bst = {}
        self.bend = {}
        for chrom in self.chrlist:
            iind = self.indf['chr'] == chrom
            bind = self.bkgdf['chr'] == chrom
            self.ist[chrom] = self.indf.loc[iind, 'start'].to_numpy()
            self.iend[chrom] = self.indf.loc[iind, 'end'].to_numpy()
            self.bst[chrom] = self.bkgdf.loc[bind, 'start'].to_numpy()
            self.bend[chrom] = self.bkgdf.loc[bind, 'end'].to_numpy()
            print(chrom, len(self.ist[chrom]), len(self.bst[chrom]))
        print("[STATUS] Finished reading inputs in " +
              str(round(time.time() - t1)) + 's')

    def read_motif_info(self):
        print("[STATUS] Reading motif set")
        self.motifs = []
        with open(self.motif_list, 'r') as f:
            for line in f:
                line = line.rstrip()
                self.motifs.append(line)

    def count_motifs(self):
        for motif in self.motifs:
            self.count_single_motif(motif)

    def count_single_motif(self, motif):
        print("[STATUS] Motif:", motif)
        mfile = self.motif_pref + motif + ".tsv.gz"
        print("[STATUS] Reading from", mfile)
        # names=['chr','start','strand', 'set','score'])
        # TODO: Streaming it sorted might be faster...
        self.ict = {}
        self.bct = {}
        t1 = time.time()
        with gzip.open(mfile, 'r') as f:
            for line in f:
                line = line.decode("utf-8")
                line = line.rstrip().split("\t")
                loc = int(line[1])
                mset = line[3]
                chrom = line[0]
                if mset not in self.ict.keys():
                    self.ict[mset] = 0
                    self.bct[mset] = 0
                self.ict[mset] = self.ict[mset] + (sum(
                    (loc<=self.iend[chrom]) * (loc>=self.ist[chrom])) > 0)
                self.bct[mset] = self.bct[mset] + (sum(
                    (loc<=self.bend[chrom]) * (loc>=self.bst[chrom])) > 0)
        print("[STATUS] Finished reading inputs in " +
              str(round(time.time() - t1)) + 's')
        print(self.ict)
        print(self.bct)


    # TODO FIX:
    def write_enrichments(self):
        self.pfmdf.to_csv(self.outfile, sep="\t", index=False, float_format="%0.4e")


if __name__ == "__main__":
    fire.Fire(count_motif_instances)
