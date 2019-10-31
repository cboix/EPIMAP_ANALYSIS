#!/usr/bin/python
# --------------------------------------------------
# Utility to calculate distance between
# a track and a set of tracks in relevant regions
#
# Modes: 1) Toggle use of a fixed set of regions
# 2) Toggle use of top N=20k regions between the two
#
# NOTE: Should not toggle both off,
# will be more efficient in another script.
#
# Also outputs the signal distribution
# in a fixed set of regions
# --------------------------------------------------
import os
# import re
import gzip
import socket
import time
import numpy as np
import six.moves.cPickle as pickle
from scipy.stats import spearmanr
from scipy.spatial import distance


# For cli usage:
import fire

rng = np.random


class dist_bws(object):
    def __init__(self, maintrack, tracksfile, output, mark,
                 nregions=20000, fixed_region=True):
        # Arguments:
        self.maintrack = maintrack
        self.tracksfile = tracksfile
        self.output = output
        self.mark = mark
        self.fixed_region = fixed_region
        self.nregions = nregions
        # Check domain and set home and directories::
        self.domain = socket.getfqdn()
        if 'broadinstitute.org' in self.domain:
            self.HOME = '/broad/compbio/cboix'
        else:
            self.HOME = os.getenv("HOME")
        self.datadir = self.HOME + '/EPIMAP_ANALYSIS/db/'
        self.trackdir = self.datadir + 'ChromImpute/imputed/'
        self.subsetdir = self.datadir + 'ChromImpute/subsetted_tracks/'
        self.indexdir = self.datadir + "/state_regionlists/"
        self.indexpref = self.indexdir + "imputed12model"
        self.modelfile = self.datadir + "Annotation/" + \
            "ChromHMM_model_imputed12_25_states.txt"
        self.bincutfile = self.datadir + "ChromImpute/" + \
            "binarization_cutoffs_match_avg_distr.tsv"
        self.start = time.time()
        # Initialize dictionaries for data, outputs, and regions:
        if self.mark == 'DNase-seq' or self.mark == 'ATAC-seq':
            self.altmark = 'DNase'
        elif self.mark == 'H2AFZ':
            self.altmark = 'H2A.Z'
        else:
            self.altmark = self.mark
        self.chrlist = ['chr' + str(i+1) for i in range(22)] + ['chrX']
        self.distance_dict = {}
        self.index = {}

    def preprocess(self):
        # Read in regions:
        self.define_fixed_regions()
        # Read in the main track data, intersect with regions:
        # NOTE: Automatically writes back up to cpfile:
        self.maindata = self.read_track(self.maintrack, againstfixed=True)
        print("[STATUS] Output has length: " + str(len(self.maindata)))
        print("Preprocessing data took: " +
              str(round(time.time() - self.start, 2)) + "s")

    def main(self):
        # Get cutoff for distances:
        self.get_binarization_cutoff()
        # Read in regions:
        if self.fixed_region:
            self.define_fixed_regions()
        # Read in the main track data:
        self.maindata = self.read_track(self.maintrack, self.fixed_region)
        self.ismainobs = (self.maintrack[0:5] == "FINAL")
        # Read in other track prefixes and compare all to main track:
        self.read_tracksfile()
        for self.vstrack in self.trackslist:
            self.isvsobs = (self.vstrack[0:5] == "FINAL")
            print("[STATUS] Evaluating against " + str(self.vstrack))
            self.vsdata = self.read_track(self.vstrack, self.fixed_region)
            self.calculate_distance_metrics(self.vstrack)
        # Write output of calculated slope:
        self.write_distances()
        print("Regression + writing took: " +
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
                clist = clist * 8
                cexplist = np.sort(np.concatenate([clist, clist + 1,
                                                   clist + 2, clist + 3,
                                                   clist + 4, clist + 5,
                                                   clist + 6, clist + 7]))
                idlist = np.concatenate([idlist, cexplist])
            # Store list as np.array (for indexing) for each chromosome:
            idlist = np.unique(np.sort(idlist))
            print("[STATUS] Index length: " + str(len(idlist)) + "\t" +
                  "Maximum index value: " + str(max(idlist)))
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

    def get_binarization_cutoff(self):
        # Get any states with >0.80 emission and top state
        print("[STATUS] Choosing relevant states for " + self.mark +
              " from " + self.bincutfile)
        self.cutoff = 2
        with open(self.bincutfile, "r") as f:
            for line in f:
                line = line.split("\n")[0].split("\t")
                if line[0] == self.mark:
                    self.cutoff = float(line[1])
        print("[STATUS] Cutoff for binary distances is: " + str(self.cutoff))

    def read_track(self, track, againstfixed):
        self.readstart = time.time()
        if againstfixed:
            try:
                cpfile = self.subsetdir + track + ".cp.gz"
                data = load_pickle(cpfile)
                print("[STATUS] Read from cpickle")
            except:
                data = self.read_from_bw(track, againstfixed)
                print("[STATUS] Saving to cpickle")
                save_pickle(cpfile, data)
        else:
            data = self.read_from_bw(track, againstfixed)
        print("[STATUS] Reading track took: " +
              str(round(time.time() - self.readstart, 2)) + "s")
        return(data)

    def read_from_bw(self, track, againstfixed):
        print("[STATUS] Reading in track " + track + " from gzipped BW")
        data = np.array([])
        for chrom in self.chrlist:
            print("[STATUS] Reading " + str(chrom))
            prefix = self.trackdir + chrom + "_"
            bwfile = prefix + track + ".wig.gz"
            if againstfixed:
                [self.comp, h1, h2] = read_bw_avg_idlist(
                    bwfile, self.index[chrom], window=8)
            else:
                [self.comp, h1, h2] = read_bw(bwfile)
            data = np.concatenate([data, self.comp])
        return(data)

    def read_tracksfile(self):
        self.trackslist = []
        with open(self.tracksfile, "r") as f:
            for line in f:
                self.trackslist.append(line.split("\n")[0])
        print("[STATUS] Will compare against " + str(len(self.trackslist)) +
              " other tracks.")

    def calculate_distance_metrics(self, vstrack):
        # Subset to top N regions if necessary:
        if self.nregions > 0:
            self.maxdata = np.maximum(self.maindata, self.vsdata)
            self.topidx = np.argpartition(self.maxdata,
                                          -self.nregions)[-self.nregions:]
            self.x = self.maindata[self.topidx]
            self.y = self.vsdata[self.topidx]
        else:
            self.x = self.maindata
            self.y = self.vsdata
        # Calculate normalized euclidean distance:
        normdist = 0.5 * (np.var(self.x - self.y) /
                          (np.var(self.x) + np.var(self.y)))
        # Pearson correlation:
        cormat = np.corrcoef(self.x, self.y)
        cor = cormat[0, 1]
        # Spearman correlation:
        rho, pval = spearmanr(self.x, self.y)
        # Also calculate binary:
        if self.ismainobs:
            self.bx = 1.0 * (self.x >= 2.0)
        else:
            self.bx = 1.0 * (self.x >= self.cutoff)
        if self.isvsobs:
            self.by = 1.0 * (self.y >= 2.0)
        else:
            self.by = 1.0 * (self.y >= self.cutoff)
        # Jaccard distance:
        jacc = distance.jaccard(self.bx, self.by)
        # Rogers-Tanimoto distance:
        rtd = distance.rogerstanimoto(self.bx, self.by)
        # Aggregate all:
        self.distance_dict[vstrack] = [cor, normdist, rho, jacc, rtd]
        print(self.distance_dict[vstrack])

    def write_distances(self):
        # Write dataframe to the region_distance directory:
        with open(self.output, "w") as f:
            for key in self.distance_dict.keys():
                vals = self.distance_dict[key]
                row = [str(self.maintrack),
                       str(key),
                       str(self.nregions),
                       str(self.fixed_region)]
                sv = [str(round(v, 4)) for v in vals]
                row = row + sv
                f.write("\t".join(row) + "\n")


def read_bw(file):
    array = []
    with gzip.open(file, 'rt') as f:
        # Read/skip two header lines:
        header1 = next(f)
        header2 = next(f)
        for line in f:
            line = line.strip("\n")
            array.append(float(line))
    narray = np.array(array)
    return([narray, header1, header2])


def read_bw_idlist(file, idlist):
    # NOTE: idlist must be sorted/unique:
    if not is_sorted(idlist):
        idlist = np.unique(np.sort(idlist))
    array = []
    with gzip.open(file, 'rt') as f:
        # Read/skip two header lines:
        header1 = next(f)
        header2 = next(f)
        ind = 0
        bufsize = 65536
        while True:
            lines = f.readlines(bufsize)
            if not lines or len(idlist) == 0:
                break
            for line in lines:
                # Add only if in list, othw. keep going:
                if ind == idlist[0]:
                    line = line.strip("\n")
                    array.append(float(line))
                    idlist = idlist[1:]
                ind = ind + 1
                # Stop if no more ids to pull
                if len(idlist) == 0:
                    break
    narray = np.array(array)
    return([narray, header1, header2])


# Return the averaged windows (200bp):
def read_bw_avg_idlist(file, idlist, window=8):
    # NOTE: idlist must be sorted/unique:
    if not is_sorted(idlist):
        idlist = np.unique(np.sort(idlist))
    if (len(idlist) % window) != 0:
        print("[ERROR] IDLIST MUST BE MULTIPLE OF WINDOWSIZE")
    array = []
    with gzip.open(file, 'rt') as f:
        # Read/skip two header lines:
        header1 = next(f)
        header2 = next(f)
        ind = 0
        widx = 0
        val = 0
        bufsize = 65536
        while True:
            lines = f.readlines(bufsize)
            if not lines or len(idlist) == 0:
                break
            for line in lines:
                # Add only if in list, othw. keep going:
                if ind == idlist[0]:
                    line = line.strip("\n")
                    val = val + float(line)
                    widx = widx + 1  # Number of values added
                    # Average every eight values:
                    if widx == window:
                        array.append(val / window)
                        val = 0
                        widx = 0
                    idlist = idlist[1:]
                ind = ind + 1
                # Stop if no more ids to pull
                if len(idlist) == 0:
                    break
    narray = np.array(array)
    return([narray, header1, header2])


# Compressed save:
def save_pickle(filename, array):
    print("[STATUS] Saving pickled to " + filename)
    try:
        with gzip.open(filename, 'wb') as outfile:
            pickle.dump(array, outfile, protocol=2)
    except:
        with gzip.open(filename, 'wb') as outfile:
            pickle.dump(array, outfile, protocol=4)


def load_pickle(filename):
    print("[STATUS] Loading pickled from " + filename)
    with gzip.open(filename, 'rb') as infile:
        array = pickle.load(infile)
    return array


# Lazy/fast check of sortedness
def is_sorted(a):
    for i in range(a.size-1):
        if a[i+1] < a[i]:
            return(False)
    return(True)


def write_bw(array, header, file):
    with gzip.open(file, 'wb') as f:
        # Write two header lines:
        f.write(header[0])
        f.write(header[1])
        for item in array:
            f.write(str(round(item, 3)) + "\n")

if __name__ == "__main__":
    fire.Fire(dist_bws)
