#!/usr/bin/python
# --------------------------
# Utility to prune lead SNPs
# --------------------------
# import os
# import re
# import socket
import gzip
import time
import numpy as np

# For cli usage:
import fire

rng = np.random


class prune_leadsnps(object):
    def __init__(self, snpfile, outfile, within=1e6, col=9):
        # Arguments:
        self.snpfile = snpfile
        self.outfile = outfile
        self.within = within
        self.col = col - 1

    def main(self):
        self.start = time.time()
        self.read_file()
        self.prune_snps()
        self.write_snps()
        print("Pruning + writing took: " +
              str(round(time.time() - self.start, 2)) + "s")

    def read_file(self):
        print("Reading snps")
        self.clist = []
        self.llist = []
        self.plist = []
        with gzip.open(self.snpfile, 'rb') as f:
            for line in f:
                line = line.strip("\n")
                line = line.split("\t")
                if line[0] == '23':
                    line[0] = 'X'
                self.clist.append(line[0])
                self.llist.append(int(line[1]))
                self.plist.append(float(line[self.col]))
        print("Read in " + str(len(self.plist)) + " snps")
        self.traversal = np.argsort(self.plist)
        self.chrlist = np.unique(self.clist)
        self.llist = np.array(self.llist)

    def prune_snps(self):
        print("Pruning snps with dist: " + str(self.within))
        self.cdict = {}
        self.pdict = {}
        # Initialize kept snps:
        j = 0
        for chrom in self.chrlist:
            self.cdict[chrom] = np.array([])
            self.pdict[chrom] = []
        # Greedily prune snps:
        for i in self.traversal:
            chrom = self.clist[i]
            if len(self.cdict[chrom]) == 0:
                dist = self.within * 10
            else:
                dist = np.min(abs(self.cdict[chrom] - self.llist[i]))
            if dist > self.within:
                j = j + 1
                self.cdict[chrom] = np.append(self.cdict[chrom], self.llist[i])
                self.pdict[chrom].append(self.plist[i])
                # Report out - chrom, loc, uid is sufficient:
        print("Kept in " + str(j) + " snps")
        # kept.snps = ldply(ut, df=gwdf, prune.snps)

    def write_snps(self):
        with gzip.open(self.outfile, 'wb') as f:
            for chrom in self.chrlist:
                for i in range(len(self.cdict[chrom])):
                    line = str(chrom) + "\t" + str(self.cdict[chrom][i]) + \
                        "\t" + str(self.pdict[chrom][i])
                    f.write(line + "\n")


if __name__ == "__main__":
    fire.Fire(prune_leadsnps)
