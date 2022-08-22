#!/usr/bin/env python
# Author: Benjamin T. James
import sys
import os
import gzip
import pickle
import numpy as np
from scipy.sparse import find
from sample_grep import for_sample_mark, get_location_table

def for_sample(sample_name, sample_index, mark_loc_table, mark_gz_list):
    print("For sample %s" % sample_name)
    if not os.path.isdir(sample_name):
        sys.exit(1)
    if not os.path.isfile(os.path.join(sample_name, "rna.txt")):
        sys.exit(1)
    for mark_gz in mark_gz_list:
        mark_name = os.path.basename(mark_gz).split('_')[0]
        with open(os.path.join(sample_name, mark_name + ".mark"), "w") as W:
            with gzip.open(mark_gz, 'r') as R:
                for_sample_mark(sample_index, mark_loc_table, R, W)

def run2(sample_file, mark_location_bed, mark_gz_list):
    mark_loc_table = get_location_table(mark_location_bed, function=int)
    sample_num_list = []
    sample_name_list = []
    with open(sample_file, 'r') as R1:
        for rline in R1:
            sample = rline.strip().split()
            sample_name_list.append(sample[0])
            sample_num_list.append(int(sample[1]))
    # For each mark file:
    for mark_gz in mark_gz_list:
        mark_name = os.path.basename(mark_gz).split('_')[0]
        print("working on", mark_name)
        sys.stdout.flush()
        #
        sample_mark_list = {s: open(os.path.join(s, mark_name + ".mark"), "w")
                            for s in sample_name_list}
        # For the mtx matrix, split into
        with gzip.open(mark_gz, 'r') as R:
            for rline in R:
                L = rline.decode('utf-8').strip().split()
                if L[0][0] == '%':
                    continue
                loc_index = int(L[0])
                sample_index = int(L[1])
                value = float(L[2])
                if sample_index in sample_num_list and loc_index in mark_loc_table:
                    loc = mark_loc_table[loc_index]
                    sample_name = sample_name_list[sample_num_list.index(sample_index)]
                    print("%s\t%d\t%d\t%f" % (loc[0], loc[1], loc[2], value),
                          file=sample_mark_list[sample_name])
        for k, v in sample_mark_list.items():
            v.close()


def run_on_cp(sample_file, mark_location_bed, markpref_list):
    # Get the masterlist (in dictionary format)
    mark_loc_table = get_location_table(mark_location_bed, function=int)
    # Load in the samples that we are running on:
    sample_num_list = []
    sample_name_list = []
    with open(sample_file, 'r') as R1:
        for rline in R1:
            sample = rline.strip().split()
            sample_name_list.append(sample[0])
            sample_num_list.append(int(sample[1])-1)
    print("Names:", sample_name_list)
    print("Indices:", sample_num_list)
    # For each mark prefix (csr and attr file pair):
    for markpref in markpref_list:
        markpref = markpref.replace("_csr.cp.gz", "")
        mark_name = os.path.basename(markpref).split('_')[0]
        print("working on", mark_name)
        sys.stdout.flush()
        # Get the names, check that they match the indices:
        mark_names = get_names(markpref)
        mark_names = ["BSS" + m.split("BSS")[1] for m in mark_names]
        mark_names = [m.split("_")[0] for m in mark_names]
        mark_names = [m.split(" ")[0] for m in mark_names]
        print(mark_names[0:10])
        # Print indices:
        mark_num_list = [mark_names.index(s) for s in sample_name_list]
        print(mark_num_list)
        print(sample_num_list)
        print(mark_num_list == sample_num_list)
        # Load data and subset indices:
        print("Loading data")
        X = get_data(markpref)
        print("Subsetting data")
        tmpX = X[:, mark_num_list]
        # For each sample, get their non-zero elements and write out file:
        for i, sample in enumerate(sample_name_list):
            print("Running", i, sample)
            # Make matrix for the sample
            Xcol = tmpX[:, i] # Sample values
            nnzX = find(Xcol)
            loc_index = nnzX[0]
            value = nnzX[2]
            print("Indexes:", loc_index[0:5])
            print("Values:", value[0:5])
            # Keep only if in the loc table (enhancers):
            kept_ind = [j for j, loc in enumerate(loc_index) if loc in mark_loc_table]
            if len(kept_ind) > 0:
                loc_index = loc_index[kept_ind]
                value = value[kept_ind]
                value = [round(v,5) for v in value]
                print("Subset to", len(kept_ind), "enhancers")
                locmat = np.array([mark_loc_table[loc] for loc in loc_index])
                # Turn into np array, bind with value:
                print(locmat.shape)
                samplemat = np.hstack((locmat, np.array(value)[:,np.newaxis]))
                print(samplemat.shape)
                print(samplemat[0:5,])
                with open(os.path.join(sample, mark_name + ".mark"), 'w') as f:
                    preformatted_write(samplemat, f, fmtstring="%s\t%s\t%s\t%s")


def run(sample_file, mark_location_bed, mark_gz_list):
    mark_loc_table = get_location_table(mark_location_bed, function=int)
    with open(sample_file, "r") as R1:
        for rline in R1:
            sample = rline.strip().split()
            sample_name = sample[0]
            sample_num = int(sample[1])
            for_sample(sample_name, sample_num, mark_loc_table, mark_gz_list)


# Faster than for loop or np.savetxt:
def preformatted_write(mat, f, fmtstring=None):
    if fmtstring is None:
        fmtstring = '\t'.join(['%g']*mat.shape[1])
    fmt = '\n'.join([fmtstring]*mat.shape[0])
    data = fmt % tuple(mat.ravel())
    f.write(data)


def get_data(prefix):
    datafile = prefix + "_csr.cp.gz"
    if not os.path.isfile(datafile):
        print(datafile + " does not exist");
        sys.exit(1)
    with gzip.open(datafile, 'rb') as f:
        X = pickle.load(f, encoding='latin1')
    return(X)


def get_names(prefix):
    attrfile = prefix + "_attr.cp.gz"
    if not os.path.isfile(attrfile):
        print(attrfile + " does not exist");
        sys.exit(1)
    with gzip.open(attrfile, 'rb') as f:
        attr = pickle.load(f, encoding='latin1')
    return(attr['names'])



if __name__ == "__main__":
    if len(sys.argv[1:]) < 3:
        print("Usage: %s sample_file mark_location_bed mark_1_gz mark_2_gz .." % sys.argv[0])
        sys.exit(1)
    # run2(sys.argv[1], sys.argv[2], sys.argv[3:])
    run_on_cp(sys.argv[1], sys.argv[2], sys.argv[3:])
    # for_sample(sys.argv[1], sys.argv[3], sys.argv[4], sys.argv[5:])
