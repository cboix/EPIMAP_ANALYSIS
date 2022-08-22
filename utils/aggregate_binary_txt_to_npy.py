#!/usr/bin/python
# -------------------------------------------------------------------
# Script to turn all txt files in directory into sparse binary files.
# Take txt.gz or wig.gz files and save them as sparse matrices
# While collating multiple sources
# -------------------------------------------------------------------
import os
import fire
import csv
import gzip
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix, hstack, vstack
# import argparse
from tqdm import tqdm

# Load config into os.environ:
import parse_config
import auxfunc_masterlist as AUX


# TODOLIST for aggregate with ML from WM:
# ---------------------------------------
# Enable only using the main 833 ids
# Add read in DHS loc (by chmm bin)
# export DML_DIR=${DBDIR}/DHS_Index_WM201902/
# Add only keep loc in DHS bins
# TODO: Add make 1 column for each state x cell
# TODO: Enable creating equivalent matrix for DNase
#       so we can take the intersection
# TODO: Second utility to perform intersection with
#       Either DNase matrix or the membership matrix

# Uses infofile:
# sample    prefix  suffix
# Allows reading in 25 or 200bp
# Allows extension off DNase-seq EXTbp (label out file)
# TODO: Processing:
# Aggregate three types of files (from masterlist):
# 1. 200bp Annotations (Binary, avgd)
# 2. 25bp H3K27ac signal (avg log10pval)
# 3. 200bp DNase-seq (Binary - keep internal) (TODO: use 25bp?)

class convert_binary(object):
    def __init__(self, infofile, dir, out, chromhmm=False, states=0,
                 intindex=True, chrom=None, onlyqc=True,
                 resolution=200, extend=0, mergestates=False, verbose=True,
                 nonovl=False, sparsify=True, idxpref=None):
        # Arguments:
        self.infofile = infofile  # Info table - what to run
        self.dir = dir  # Directory where files located
        self.out = out  # Output prefix
        self.chromhmm = chromhmm
        self.states = states
        if self.chromhmm:
            print("STATES: " + str(self.states))
        self.intindex = intindex  # Intersect w/ index:
        self.chrom = chrom
        self.onlyqc = onlyqc  # Keep only samples passing our qc
        self.resolution = int(resolution)
        self.extend = int(extend)
        self.mergestates = mergestates
        self.nonovl = nonovl
        self.sparsify = sparsify  # Controls if we make sparse or dense matrices
        if self.resolution == 200:
            self.sparsify = True  # Always bin the 200 bp res tracks
        self.verbose = verbose
        if self.chrom is None:
            self.chrlist = ['chr' + str(i+1) for i in range(22)] + ['chrX']
        else:
            if type(chrom) == str:
                self.chrlist = [chrom]
            else:
                self.chrlist = chrom
        # Files and directories:
        self.cpdir = self.dir.rstrip("/") + "/cpfiles/"
        if not os.path.exists(self.cpdir):
            os.mkdir(self.cpdir)
        self.datadir = os.environ['DBDIR']
        # For files, flag with resolution/extension:
        self.re_mid = "_r" + str(self.resolution) + "_e" + str(self.extend)
        self.keptlist = os.environ['ANNDIR'] + '/kept_bssid_20190322.txt'
        # Masterlist locations:
        if idxpref is not None:
            self.dmlpref = idxpref
            self.flag = '_alt'
        else:
            if self.nonovl:
                print("[STATUS] Using the non-overlapping DHS list:")
                self.dmlpref = os.environ['NON_DMLPREF']
                print(self.dmlpref)
                self.flag = '_nonovl'
            else:
                self.dmlpref = os.environ['DMLPREF']
                self.flag = None
        self.indexlocfile = self.dmlpref + '_hg19.core.srt.txt'
        self.cplocfile = self.dmlpref + '_hg19' + \
            self.re_mid + ".core.srt.cp.gz"
        self.cpnamfile = self.dmlpref + '_hg19' + \
            self.re_mid + "_names.core.srt.tsv"
        # https://stackoverflow.com/questions/5980042/
        if self.verbose:
            def verboseprint(self, *args):
                # Print each argument separately so caller doesn't need to
                # stuff everything to be printed into a single string
                for arg in args:
                    print(arg)
                print
        else:
            self.verboseprint = lambda *a: None  # do-nothing function
        self.ord_il = False  # Order the index by bin # (False = as is)

    def main(self):
        self.setup()
        self.process()

    def setup(self):
        self.load_infofile()
        # Load index if necessary:
        if self.intindex:
            self.load_index()
        # Load samples:
        if self.onlyqc:
            self.load_keptsamples()
        # Filter infofile against sample list:
        if self.onlyqc:
            self.pid = [i for i, p in enumerate(self.ids) if p in self.samples]
        else:
            self.pid = list(range(len(self.ids)))
        print("[STATUS] There are " + str(len(self.pid)) +
              " prefixes of files in sample list.")

    def load_infofile(self):
        self.ids = []
        self.prefixes = []
        self.suffixes = []
        self.dirs = []
        with open(self.infofile, 'r') as f:
            for line in f:
                row = line.rstrip("\n").split("\t")
                self.ids.append(row[0])
                self.prefixes.append(row[1])
                self.suffixes.append(row[2])
                if len(row) > 3:
                    self.dirs.append(row[3])
                else:
                    self.dirs.append(self.dir)
        # Make any relevant cpdirs
        for dir in np.unique(self.dirs):
            cpdir = dir.rstrip("/") + "/cpfiles/"
            if not os.path.exists(cpdir):
                os.mkdir(cpdir)
        print("[STATUS] Read in infofile, with " + str(len(self.ids)) + " ids.")

    def load_index(self):
        try:
            # Try to load processed locations (faster)
            readfile = self.cplocfile
            if os.path.isfile(readfile):
                print("[STATUS] Loading processed masterlist: " + readfile)
                out = AUX.load_pickle_gzip(readfile)
                self.idict = out['loc']
                self.ndict = out['names']
                self.wdict = out['weights']
                self.ddict = out['indices']
            else:
                raise Exception("No such file: " + readfile)
        except:
            # Otherwise read the file and process it:
            print("[STATUS] Processing masterlist")
            self.ndict = {}  # Chunk names
            self.idict = {}  # Bins
            self.wdict = {}  # Weight for each bin (amt inside)
            self.ddict = {}  # Count
            self.oldchr = ''
            with open(self.indexlocfile, 'r') as f:
                count = 0
                for line in f:
                    self.process_index_line(line, count)
                    count = count + 1  # Track line number for sorting later.
            # Process the aggregate:
            for key in self.idict.keys():
                print(key)
                # Flatten:
                il = np.array(flatten_list(self.idict[key]))
                wl = np.array(flatten_list(self.wdict[key]))
                nl = np.array(flatten_list(self.ndict[key]))
                dl = np.array(flatten_list(self.ddict[key]))
                # print('[STATUS] Inds:', dl[0:20])
                # print('[STATUS] Nams:', nl[0:20])
                # Order lists according to bin #:
                if self.ord_il:
                    order = np.argsort(il)  # Original sorted order
                else:
                    # print("Ordering by num")
                    order = np.argsort(dl, kind='mergesort')  # New sorted order (as read in)
                # print('[STATUS] Order looks like:', order[0:20])
                # print('[STATUS] Ordered:', dl[order][0:20])
                self.idict[key] = il[order]
                self.wdict[key] = wl[order]
                self.ndict[key] = nl[order]
                self.ddict[key] = dl[order]
            # Save it as pickled file
            out = {'names': self.ndict, 'loc': self.idict,
                   'weights': self.wdict, 'indices': self.ddict}
            AUX.save_pickle_gzip(self.cplocfile, out)
            final_names = []
            final_chr = []
            print("[STATUS] Getting order of names:")
            for chrom in self.chrlist:
                # NOTE: np.unique does an internal sort that was causing issues
                # We must retain idx if we want to keep the order
                if self.ord_il:
                    # NOTE: Copy exactly the name sort from merge bins:
                    # TODO: Not sustainable, may cause bugs. Change later?
                    nlist = list(np.unique(self.ndict[chrom]))
                else:
                    # Order-perserving unique:
                    u, ind = np.unique(self.ndict[chrom], return_index=True)
                    nlist = list(u[np.argsort(ind, kind='mergesort')])
                final_names = final_names + nlist
                final_chr = final_chr + [chrom] * len(nlist)
            # Write out final names as readable TSV:
            print("[STATUS] Writing names list order to: " + self.cpnamfile)
            ndf = pd.DataFrame({'name': final_names,
                                'chr': final_chr,
                                'cls': list(np.arange(len(final_names)) + 1)})
            ndf.to_csv(self.cpnamfile, sep='\t', index=False)

    def process_index_line(self, line, count):
        line = line.split("\n")[0].split("\t")
        # Extend by self.extend bp in each direction:
        if self.extend != 0:
            mid = (float(line[2]) + float(line[1])) / 2.0
            n1 = mid - self.extend
            n2 = mid + self.extend
        else:
            n1 = float(line[1])
            n2 = float(line[2])
        # ChromHMM (200bp) bins with weights:
        # Full bins (round bottom up and top down)
        bin1 = int(n1 / self.resolution) + 1
        bin2 = int(n2 / self.resolution)
        fullbins = list(np.arange(bin1, bin2))
        # Fractional bins:
        frac1 = round(float(bin1) - (n1 / (1.0 * self.resolution)), 3)
        frac2 = round((n2 / (1.0 * self.resolution)) - float(bin2), 3)
        bins = [bin1 - 1] + fullbins + [bin2]
        weights = [frac1] + [1.0] * len(fullbins) + [frac2]
        # Add chr if not seen yet:
        if line[0] != self.oldchr:
            if line[0] not in self.idict.keys():
                print("[STATUS] Starting chromosome " + line[0])
                self.idict[line[0]] = []
                self.wdict[line[0]] = []
                self.ndict[line[0]] = []
                self.ddict[line[0]] = []
        # Locations, weights, and names:
        self.idict[line[0]].append(bins)
        self.wdict[line[0]].append(weights)
        self.ndict[line[0]].append([line[3]] * len(bins))
        self.ddict[line[0]].append([count] * len(bins))
        self.oldchr = line[0]

    def load_keptsamples(self):
        self.samples = []
        with open(self.keptlist, 'r') as f:
            for line in f:
                self.samples.append(line.split("\n")[0])

    # TODO: Also process leaving ALL STATES in the matrix (833 x 3.6M is small)
    # -- Make different utility??
    # -- Then We can filter states we care about out of it.
    def process(self):
        print("[STATUS] Processing files")
        # For each chromosome, join all files from the noted prefixes
        # Concatenate each chromosome to make main file:
        if self.chrom is None:
            self.concatenate_chrom()
        else:
            for chrom in self.chrlist:
                self.process_chrom(chrom)

    def process_chrom(self, chrom):
        print("[STATUS] Processing " + chrom)
        chrompref = self.out + "_" + chrom + self.re_mid
        if self.mergestates:
            chrompref = chrompref + "_merged"
        chrdata_file = chrompref + "_csr.cp.gz"
        chrattr_file = chrompref + "_attr.cp.gz"
        j = 0
        namelist = []
        for idnum in tqdm(self.pid):
            currid = self.ids[idnum]
            filepref = self.prefixes[idnum] + chrom + self.suffixes[idnum]
            maindir = self.dirs[idnum]
            self.verboseprint(currid + ": " + filepref)
            if self.intindex:
                # Get intersected matrix:
                [X, names] = self.get_mat_idlist(filepref, chrom, maindir)
            else:
                # Get full matrix for states
                [X, names] = self.get_mat(filepref, maindir)
            if type(names) != list:
                names = [names]
            # Concatenate:
            if j == 0:
                FULL = X
                namelist = names
                j = 1
            else:
                if self.sparsify:
                    FULL = hstack([FULL, X])
                else:
                    FULL = np.hstack([FULL, X])
                namelist = namelist + names
                # print("[STATUS] Current shape: " + str(FULL.shape))
                # print("[STATUS] Current names: " + str(namelist))
        # Print out dataset:
        print("[STATUS] Writing chromosome " + chrom)
        self.save_cpfiles(FULL, namelist, main=chrdata_file, attr=chrattr_file)

    # Standalone to preprocess files:
    def preprocess(self, prefix):
        self.setup()
        idnumlist = np.where(np.array(self.ids) == prefix)[0]
        for idnum in idnumlist:
            print(idnum)
            currid = self.ids[idnum]
            maindir = self.dirs[idnum]
            print("[STATUS] Preprocessing row " + str(idnum) +
                  " of samplesheet.")
            for chrom in self.chrlist:
                print
                print("[STATUS] Processing " + prefix +
                    " (" + currid + ") for " + chrom)
                filepref = self.prefixes[idnum] + chrom + self.suffixes[idnum]
                self.verboseprint(currid + ": " + filepref)
                # Intersect with DHS index + current run
                # Saves both intersect and further files for later loading
                [X, names] = self.get_mat_idlist(filepref, chrom, maindir)

    # Get the raw intersection of bins vs. txt.gz file:
    # NOTE: Works for both dense and sparse matrix
    def get_intersected_file(self, filepref, chrom, maindir):
        print("[STATUS] Getting intersection of file and DHS list")
        idlist = list(self.idict[chrom])  # Is a sorted IDLIST
        fullname = maindir + filepref
        cppref = maindir + "/cpfiles/" + filepref.split(".gz")[0] + self.re_mid
        if self.flag is not None:
            cppref = cppref + self.flag
        cp_file = cppref + '_intersect_csr.cp.gz'
        attr_file = cppref + '_intersect_attr.cp.gz'
        try:
            [X, names] = self.load_cpfiles(main=cp_file, attr=attr_file)
        except:
            print("Reading: " + fullname)
            with gzip.open(fullname, 'rt') as f:
                header = next(f)
                header2 = next(f)
                names = header.split("\t")[0]
                X = []
                ind = 0
                for line in f:
                    X.append(line)
                    ind += 1
                if self.resolution == 200:
                    tmpX = [int(X[i - 1].split("\n")[0]) for i in idlist]
                else:
                    tmpX = [float(X[i - 1].split("\n")[0]) for i in idlist]
                X = np.array(tmpX)
            self.save_cpfiles(X, names, main=cp_file, attr=attr_file)
        print("[STATUS] Extracted " + str(len(X)) + " bins. There are " +
              str(len(self.idict[chrom])) + " index bins.")
        return([X, names])

    # Get matrix in the context of an index list:
    def get_mat_idlist(self, filepref, chrom, maindir):
        # Define CP file names by task:
        cppref = maindir + "/cpfiles/" + filepref.split(".gz")[0] + self.re_mid
        if not self.sparsify and self.resolution == 25:
            cppref = cppref + '_dense'
        # TODO: Fix nonovl to be more general descriptor of list:
        if self.flag is not None:
            cppref = cppref + self.flag
        if self.chromhmm:
            if self.states == 0:
                cppref = cppref + "_dhsidx_allstates"
            else:
                st = [str(state) for state in self.states]
                cppref = cppref + "_dhsidx_" + "_".join(st)
            if self.mergestates:
                cppref = cppref + "_merged"
        cp_file = cppref + '_csr.cp.gz'
        attr_file = cppref + '_attr.cp.gz'
        self.verboseprint(cp_file)
        # Load files:
        try:
            [X, names] = self.load_cpfiles(main=cp_file, attr=attr_file)
        except:
            [X, names] = self.get_intersected_file(filepref, chrom, maindir)
            if self.resolution == 25:
                [X, chunks] = self.merge_bins_average(X, chrom)
            else:
                # Filter to 0/1 if states:
                X = self.binarize_matrix(X, tocsr=False)
                # Merge into bins:
                [X, chunks] = self.merge_bins(X, chrom)
                # NOTE:
            if self.chromhmm and not self.mergestates:
                names = [names + "_" + str(s) for s in self.states]
            print("[STATUS] Merged into " + str(len(chunks)) + " chunks.")
            self.save_cpfiles(X, names,
                              main=cp_file, attr=attr_file, chunks=chunks)
        return([X, names])

    # Merge the bins based on frac, maximal count:
    # Currently handle ties with taking first state
    # TODO: Merge the bins 0/1 first.
    def merge_bins(self, X, chrom):
        print("[STATUS] Merging bins:")
        # NOTE: X can be either 0/1 or 1-25 ish
        # Sort by names:
        # TODO: This name sort has to be carried on later!
        if self.ord_il:
            order = np.argsort(self.ndict[chrom])
        else:
            order = np.argsort(self.ddict[chrom], kind='mergesort')
        wlist = self.wdict[chrom][order]
        nlist = self.ndict[chrom][order]
        # NOTE: X multicol:
        if self.chromhmm and not self.mergestates:
            Xsort = X[order, :]
        else:
            Xsort = X[order, np.newaxis]
        # Init outputs:
        chunks = []
        chunk = nlist[0]
        # for each column:
        NCOLS = Xsort.shape[1]
        st = [0] * NCOLS
        NNAM = len(np.unique(nlist))
        out = np.zeros((NNAM, NCOLS), dtype=np.int8)
        if self.chromhmm and self.states == 0:
            arrshape = (NCOLS, 26)
        else:
            arrshape = (NCOLS, 2)
        arr = np.zeros(arrshape)
        # Go through in name order:
        ni = 0
        NLEN = len(nlist)
        for i in range(NLEN):
            n = nlist[i]
            # Merge chunks:
            if n != chunk or i == (NLEN - 1):
                for j in range(NCOLS):
                    a = np.argmax(arr[j, :])
                    # NOTE: Works best for mergestates
                    # Only if majority fraction:
                    if arr[j,a] > 0.5:
                        st[j] = a
                    else:
                        st[j] = 0
                # out.append(st)
                chunks.append(chunk)
                out[ni, :] = st
                # Reset values:
                chunk = n
                ni = ni + 1
                arr = np.zeros(arrshape)
            # Add new values:
            for j in range(NCOLS):
                v = Xsort[i, j]
                w = wlist[i]
                arr[j, v] = arr[j, v] + w
        out = np.array(out)
        print("NNZ: " + str(np.sum(out, axis=0)) +
              " out of " + str(out.shape[0]))
        # NOTE: Binary should always be outputted as sparse:
        csarr = csr_matrix(out)
        return([csarr, chunks])

    # Merge the bins by weighted average
    def merge_bins_average(self, X, chrom):
        print("[STATUS] Merging bins:")
        # NOTE: X can be either 0/1 or 1-25 ish
        # Sort by names:
        if self.ord_il:
            order = np.argsort(self.ndict[chrom])
        else:
            order = np.argsort(self.ddict[chrom], kind='mergesort')
        wlist = self.wdict[chrom][order]
        nlist = self.ndict[chrom][order]
        # NOTE: X multicol:
        order = np.array(order)
        X = np.array(X)
        if self.chromhmm and not self.mergestates:
            Xsort = X[order, :]
        else:
            Xsort = X[order, np.newaxis]
        print(Xsort.shape)
        # Init outputs:
        chunks = []
        chunk = nlist[0]
        # for each column:
        NCOLS = Xsort.shape[1]
        st = [0.0] * NCOLS
        NNAM = len(np.unique(nlist))
        out = np.zeros((NNAM, NCOLS), dtype=np.float64)
        # Go through in name order:
        ni = 0
        NLEN = len(nlist)
        val = 0.0
        weight = 0.0
        for i in range(NLEN):
            n = nlist[i]
            # Merge chunks:
            if n != chunk or i == (NLEN - 1):
                for j in range(NCOLS):
                    # TODO: Indexed by j
                    # Averaged value:
                    a = float(val / weight)
                    st[j] = a
                chunks.append(chunk)
                out[ni, :] = st
                # Reset values:
                chunk = n
                ni = ni + 1
                val = 0.0
                weight = 0.0
            # Add new values:
            for j in range(NCOLS):
                v = Xsort[i, j]
                w = wlist[i]
                val = val + v
                weight = weight + w
        out = np.array(out)
        print("Out shape: " + str(out.shape[0]))
        # NOTE: Sparsify by cutting out values less than 2:
        if self.sparsify:
            out[out < 2.0] = 0
            csarr = csr_matrix(out)
            return([csarr, chunks])
        else:
            # Return the dense matrix otherwise:
            return([out, chunks])

    # Get matrix agnostic of any index list:
    def get_mat(self, filepref, maindir):
        # Names of files raw/cp:
        fullname = maindir + filepref
        cppref = maindir + "/cpfiles/" + filepref.split(".gz")[0]
        if self.flag is not None:
            cppref = cppref + self.flag
        if self.chromhmm:
            st = [str(state) for state in self.states]
            cppref = cppref + "_" + "_".join(st)
        cp_file = cppref + '_csr.cp.gz'
        attr_file = cppref + '_attr.cp.gz'
        # Load files
        try:
            [X, names] = self.load_cpfiles(main=cp_file, attr=attr_file)
        except:
            with gzip.open(fullname, 'rt') as handle:
                reader = csv.reader(handle, delimiter='\t')
                header = next(reader)
                names = next(reader)
                X = []
                for line in reader:
                    row = [float(x) for x in line]
                    X.append(row)
            X = self.binarize_matrix(X)
            # Save final matrix:
            self.save_cpfiles(X, names, main=cp_file, attr=attr_file)
        return([X, names])

    def binarize_matrix(self, X, tocsr=True):
        X = np.array(X)
        if self.chromhmm and self.states != 0:
            if self.mergestates:
                print("Merging states")
                ixmat = 1 * np.in1d(X, np.array(self.states))
                print(ixmat.shape)
                print("NNZ: " + str(np.sum(ixmat)) +
                      " out of " + str(len(ixmat)))
            else:
                ixmat = 1 * np.array([X == state for state in self.states]).T
                print(ixmat.shape)
                print("NNZ: " + str(np.sum(ixmat, axis=0)) +
                      " out of " + str(ixmat.shape[0]))
            if tocsr:
                mat = csr_matrix(ixmat)
            else:
                mat = np.array(ixmat)
        else:
            if tocsr:
                mat = csr_matrix(X)  # NOTE: May not be efficient for states=0
            else:
                mat = np.array(X)
        return(mat)

    def save_cpfiles(self, X, names, main, attr, chunks=None):
        print("[STATUS] Saving to: " + main)
        out = {'names': names}
        if chunks is not None:
            out['chunks'] = chunks
        AUX.save_pickle_gzip(main, X)
        AUX.save_pickle_gzip(attr, out)

    def load_cpfiles(self, main, attr):
        if os.path.isfile(main):
            self.verboseprint("[STATUS] Trying to load cp: " + main)
            X = AUX.load_pickle_gzip(main)
            out = AUX.load_pickle_gzip(attr)
            names = out['names']
            self.verboseprint("[STATUS] Loaded cp file successfully")
        else:
            raise Exception("No such file: " + main)
        return([X, names])

    def concatenate_chrom(self):
        print("[STATUS] Concatenating all chromosomes")
        out_pref = self.out + self.re_mid + "_allchr"
        if self.mergestates:
            out_pref = out_pref + "_merged"
        out_attr = out_pref + "_attr.cp.gz"
        full_names = None
        self.X = None
        for chrom in tqdm(self.chrlist):
            chrompref = self.out + "_" + chrom + self.re_mid
            if self.mergestates:
                chrompref = chrompref + "_merged"
            chrdata_file = chrompref + "_csr.cp.gz"
            chrattr_file = chrompref + "_attr.cp.gz"
            # data_file = self.out + "_" + chrom + self.re_mid
            # attr_file = data_file + "_attr.cp.gz"
            if not os.path.isfile(chrdata_file):
                self.process_chrom(chrom)
            if self.sparsify:
                X_chr = csr_matrix(AUX.load_file_save_sparse(chrompref))
            else:
                X_chr = AUX.load_file_save_sparse(chrompref)
            attr = AUX.load_pickle_gzip(chrattr_file)
            if full_names is None:
                full_names = attr['names']
            if attr['names'] != full_names:
                raise ValueError("Not the same names!")
            if self.X is None:
                self.X = X_chr
            else:
                if self.sparsify:
                    self.X = vstack([self.X, X_chr])
                else:
                    self.X = np.vstack([self.X, X_chr])
        print("[STATUS] Writing final matrix")
        # Add indexes dictionary and save attributes:
        AUX.save_pickle_gzip(out_attr, attr)
        # Save as CSR sparse (both NPZ and CP)
        AUX.save_pickle_gzip(out_pref + "_csr.cp.gz", self.X)
        # Write out the name order:
        AUX.savelist(attr['names'], fname=out_pref + "_names.tsv")
        if self.sparsify:
            AUX.save_sparse_csr(out_pref + "_csr", self.X)
        else:
            # Save hdf5 matrix??
            pass
        if self.chromhmm and type(self.states) == list:
            self.write_collapsed_states_matrix(out_pref)

    # Write matrix method specific to enhancer states:
    def write_collapsed_states_matrix(self, pref):
        collfile = pref + "_collapsed_csr.cp.gz"
        NS = len(self.states)
        if not os.path.exists(collfile):
            # To COO, reduced representation:
            self.Xcoo = coo_matrix(self.X)
            col = self.Xcoo.col / NS
            col = col.astype('int')
            Xnew = coo_matrix((self.Xcoo.data, (self.Xcoo.row, col)),
                              (self.X.shape[0], int(self.X.shape[1] / NS)))
            Xcsr = csr_matrix(Xnew)
            print("[STATUS] Collapsed matrix to size: " + str(Xcsr.shape))
            AUX.save_pickle_gzip(collfile, Xcsr)
            print("[STATUS] Done writing out collapsed CSR")


def flatten_list(l):
    fl = [item for sublist in l for item in sublist]
    return(fl)


if __name__ == "__main__":
    fire.Fire(convert_binary)
