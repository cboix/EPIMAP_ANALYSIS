#!/usr/bin/python
# ------------------------------------------
# Utility to extract the signal distribution
# in a fixed set of relevant regions to mark
# ------------------------------------------
import os
# import re
import gzip
import socket
import time
import numpy as np

# For cli usage:
import fire

rng = np.random


class extract_signal_distribution(object):
    def __init__(self, maintrack, output, mark, trackdir=None):
        # Arguments:
        self.maintrack = maintrack
        self.output = output
        self.mark = mark
        # Check domain and set home and directories::
        self.domain = socket.getfqdn()
        if 'broadinstitute.org' in self.domain:
            self.HOME = '/broad/compbio/cboix'
        else:
            self.HOME = os.getenv("HOME")
        self.datadir = self.HOME + '/EPIMAP_ANALYSIS/db/'
        if trackdir is None:
            self.trackdir = self.datadir + 'ChromImpute/imputed/'
        else:
            self.trackdir = trackdir
        self.indexdir = self.datadir + "/state_regionlists/"
        self.indexpref = self.indexdir + "imputed12model"
        # TODO: Should add parameter for model:
        self.modelfile = self.datadir + "Annotation/" + \
            "ChromHMM_model_imputed12_25_states.txt"
        self.start = time.time()
        # Initialize dictionaries for data, outputs, and regions:
        if self.mark == 'DNase-seq' or self.mark == 'ATAC-seq':
            self.altmark = 'DNase'
        elif self.mark == 'H2AFZ':
            self.altmark = 'H2A.Z'
        else:
            self.altmark = self.mark
        self.chrlist = ['chr' + str(i+1) for i in range(22)] + ['chrX']
        self.distr_dict = {}
        self.state_distr_dict = {}
        self.index = {}

    def main(self):
        # Read in regions:
        self.define_fixed_regions()
        # Read in the track data and add to dict:
        self.read_track(self.maintrack)
        # Write outputs as one file (dataframe):
        self.write_distr_dataframe()
        print("Reading + writing took: " +
              str(round(time.time() - self.start, 2)) + "s")

    def define_fixed_regions(self):
        # Using mark, choose an element set from ChromHMM and pick regions
        print("[STATUS] Reading in fixed regions for mark " + str(self.mark))
        # Read in model emissions and get the appropriate state(s) for the mark
        self.choose_relevant_states()
        # Read in the predetermined fixed regions (all regions with label)
        for chrom in self.chrlist:
            idlist = np.array([])
            print("[STATUS] Reading in fixed regions for " + chrom)
            for state in self.relevant_states:
                filename = self.indexpref + "_s" + str(state) + \
                    "_" + chrom + "_possible_regions.tsv"
                clist = np.loadtxt(filename, dtype=int)
                clist = np.array(clist) * 8  # scale up by 200 / 25
                cexplist = np.sort(np.concatenate([clist, clist + 1,
                                                   clist + 2, clist + 3,
                                                   clist + 4, clist + 5,
                                                   clist + 6, clist + 7]))
                idlist = np.concatenate([idlist, cexplist])
            # Store list as np.array (for indexing) for each chromosome:
            idlist = np.unique(np.sort(idlist))
            print("[STATUS] Index length: " + str(len(idlist)))
            self.index[chrom] = idlist.astype(int)

    def choose_relevant_states(self):
        # Get any states with >0.80 emission and top state
        print("[STATUS] Choosing relevant states for " + self.mark +
              " from " + self.modelfile)
        states80 = []
        topstate = -1
        prob = 0
        with open(self.modelfile, "r") as f:
            for line in f:
                line = line.split("\n")[0].split("\t")
                if line[0] == 'emissionprobs' and line[3] == self.altmark:
                    if line[4] == '1':
                        state = line[1]
                        newprob = float(line[5])
                        if newprob > prob:
                            prob = newprob
                            topstate = state
                        # Also write to 0.80 list (topstates get replaced)
                        if newprob > 0.80:
                            states80.append(state)
        self.relevant_states = list(np.unique(np.sort(states80 + [topstate])))
        print(self.relevant_states)

    def read_track(self, track):
        self.readstart = time.time()
        print("[STATUS] Reading in track " + track)
        for chrom in self.chrlist:
            print("[STATUS] Reading " + str(chrom))
            prefix = self.trackdir + chrom + "_"
            bwfile = prefix + track + ".wig.gz"
            # NOTE: idlist must be sorted/unique:
            idlist = self.index[chrom]
            print("[STATUS] Index length: " + str(len(idlist)) + "\t" +
                  "Maximum index value: " + str(max(idlist)))
            statedata = []
            otherdata = []
            with gzip.open(bwfile, 'rt') as f:
                # Read/skip two header lines:
                line = next(f)
                line = next(f)
                ind = 0
                bufsize = 65536
                while True:
                    lines = f.readlines(bufsize)
                    if not lines:
                        break
                    for line in lines:
                        # Add to either list:
                        value = round(float(line.strip("\n")), 1)
                        if len(idlist) > 0 and ind == idlist[0]:
                            statedata.append(value)
                            idlist = idlist[1:]
                        else:
                            otherdata.append(value)
                        ind = ind + 1

            # Process the collected data:
            u, c = np.unique(statedata, return_counts=True)
            sdchr = (dict(zip(u, c)))
            self.state_distr_dict = merge_dict(self.state_distr_dict, sdchr)
            # For other loc too:
            u, c = np.unique(otherdata, return_counts=True)
            odchr = (dict(zip(u, c)))
            self.distr_dict = merge_dict(self.distr_dict, odchr)
        print("[STATUS] Reading track took: " +
              str(round(time.time() - self.readstart, 2)) + "s")

    def write_distr_dataframe(self):
        # Write distances to the region_distance directory:
        with open(self.output, "w") as f:
            keys = np.sort(self.distr_dict.keys())
            skeys = np.sort(self.state_distr_dict.keys())
            for key in keys:
                f.write(str(key) + "\t" +
                        str(round(self.distr_dict[key])) + "\t" +
                        "outside" + "\n")
            for key in skeys:
                f.write(str(key) + "\t" +
                        str(round(self.state_distr_dict[key])) + "\t" +
                        "instate" + "\n")


# https://stackoverflow.com/questions/10461531/
def merge_dict(x, y):
    merge = {k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y)}
    return(merge)


if __name__ == "__main__":
    fire.Fire(extract_signal_distribution)
