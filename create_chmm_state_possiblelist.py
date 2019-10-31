#!/usr/bin/python
# -----------------------------------------------
# Create fixed lists of possible locations of
# each chromhmm state in the imputed12 model
# ----------------------------------------------
import os
import re
import glob
import gzip
import socket
import time
import numpy as np

# For cli usage:
import fire

rng = np.random


# MUST READ IN 127 x 23 = 2921 STATEBYLINE from (or other model):
# directory='/broad/compbio/jernst/compbio-hp/' + \
#     'SIGNALPREDICT6/TIER2_REORDER_MODEL/states25/STATEBYLINE/'
# outprefix='/broad/compbio/cboix/EPIMAP_ANALYSIS/db/possible_regions/n25/imputed12model'
class create_regionlist(object):
    def __init__(self, directory, outprefix, nstate=25):
        # Arguments:
        self.directory = directory
        self.outprefix = outprefix
        self.nstate = str(nstate)
        self.suffix = "_statebyline.txt.gz"
        # Check domain and set home and directories::
        self.domain = socket.getfqdn()
        if 'broadinstitute.org' in self.domain:
            self.HOME = '/broad/compbio/cboix'
        else:
            self.HOME = os.getenv("HOME")
        self.datadir = self.HOME + '/EPIMAP_ANALYSIS/db/'
        self.chrlist = ['chr' + str(i+1) for i in range(22)] + ['chrX']
        self.start = time.time()

    def main(self):
        for chrom in self.chrlist:
            print("[STATUS] Reading/writing for " + chrom)
            self.cstart = time.time()
            self.initialize_index(chrom)
            self.list_files(chrom)
            self.read_files(chrom)
            self.write_index(chrom)
            # Write out each state + each chromosome to outdir:
            print("[STATUS] Reading/writing " + chrom + " took: " +
                  str(round(time.time() - self.cstart, 2)) + "s")
        print("[STATUS] Total time to go through all chromosomes: " +
              str(round(time.time() - self.start, 2)) + "s")

    def initialize_index(self, chrom):
        self.index = {}
        for state in range(1, int(self.nstate)):
            self.index[str(state)] = []

    def list_files(self, chrom):
        # Make list of directory:
        self.pattern = "*_" + chrom + self.suffix
        self.searchexpr = self.directory + self.pattern
        print(self.searchexpr)
        self.filelist = glob.glob(self.searchexpr)
        self.basenames = [re.sub(".*/", "", filename)
                          for filename in self.filelist]
        print("[STATUS] Adding " + str(len(self.filelist)) + " files to list.")

    def read_files(self, chrom):
        for filename in self.filelist:
            print("[STATUS] Reading " + filename)
            with gzip.open(filename, 'rt') as f:
                # Read/skip two header lines:
                line = next(f)
                line = next(f)
                ind = 0
                for line in f:
                    state = line.strip("\n")
                    if state != self.nstate:
                        self.index[state].append(ind)
                    ind = ind + 1

            # Get sorted and unique index elements:
            for state in range(1, int(self.nstate)):
                tmplist = np.unique(np.sort(self.index[str(state)]))
                self.index[str(state)] = list(tmplist)

    def write_index(self, chrom):
        for state in range(1, int(self.nstate)):
            filename = self.outprefix + "_s" + str(state) + \
                "_" + chrom + "_possible_regions.tsv"
            with open(filename, 'w') as f:
                for item in self.index[str(state)]:
                    f.write(str(item) + "\n")


if __name__ == "__main__":
    fire.Fire(create_regionlist)
