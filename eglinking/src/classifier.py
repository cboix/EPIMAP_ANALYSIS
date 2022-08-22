#!/usr/bin/env python3
# Author: Benjamin T. James
# first sample is enhancer state unless [which] is specified
import sys
import os
import time
from data import Data
from sklearn.linear_model import LogisticRegression
from scipy.stats import pearsonr
import numpy as np
import json
import pickle
import fire
from time import gmtime, strftime
import gzip
REGULARIZATION=1.0
PENALTY="l2"

# Function to make a dataset from a list of samples
def get_data(data_dir, sample_dir, sample_list, saved=None):
    # Load from pickle if saved file:
    if saved:
        print("[STATUS] Loading data from saved file:", saved)
        with gzip.open(saved, 'rb') as R:
            data = pickle.load(R)
            return data
    if not os.path.isdir(sample_dir):
        print("%s is not a valid directory" % sample_dir, file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(data_dir):
        print("%s is not a valid directory" % data_dir, file=sys.stderr)
        sys.exit(1)
    mark_list = []
    # Put together data from samples in list:
    for sample in sample_list:
        # For each sample, add the sample's marks data:
        mark_list += [f for f in os.listdir(os.path.join(sample_dir, sample))
                      if len(f) > 5 and f[-5:] == ".mark"]
    mark_list = set(mark_list)
    print("mark list:", mark_list, file=sys.stderr)
    # Create data object from all the mark datafiles:
    data = Data(data_dir, sample_dir, sample_list, mark_list)
    fname = strftime("%Y_%m_%d_%H:%M:%S", gmtime())
    with gzip.open('%s.pickle.gz' % fname, 'wb') as W:
        pickle.dump(data, W)
    return data


def kullback_leibler(A, B):
    total = 0
    for a, b in zip(A, B):
        if b != 0 and a != 0: # if a=0 then total is unchanged
            lg = np.log(a/b)
            total += a * lg
    return total

def pearson_mat(A, B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:,None]
    B_mB = B - B.mean(1)[:,None]
    # Sum of squares across rows
    ssA = (A_mA**2).sum(1);
    ssB = (B_mB**2).sum(1);
    out = (A_mA * B_mB).sum(1) / np.sqrt(ssA * ssB)
    return(out)

# Sub-function for calculating the correlation/divergence:
def calc_metric(x, expr_vals, random_expr_vals, metric):
    # for all marks, correlate sample signal vs same-sample expr
    if np.sum(x) == 0:
        pos_corr = 0
        neg_corr = 0
    else:
        if metric == 'kl':
            pos_corr = kullback_leibler(x, expr_vals)
            neg_corr = kullback_leibler(x, random_expr_vals)
        elif metric == 'pearson':
            # Catch issues with constants:
            try:
                pos_corr, pos_pval = pearsonr(x, expr_vals)
            except:
                pos_corr = 0.0
                pass
            try:
                neg_corr, neg_pval = pearsonr(x, random_expr_vals)
            except:
                neg_corr = 0.0
                pass
            if np.isnan(pos_corr):
                pos_corr = 0
            if np.isnan(neg_corr):
                neg_corr = 0
        else:
            print("[WARNING] Metric:", metric, "is not recognized/implemented")
    # print([pos_corr, neg_corr])
    return pos_corr, neg_corr


def run(data_dir=None, sample_dir=None, sample_list=None, enh_state_list=None,
        threshold=None, which_first=None, output=sys.stdout,
        saved=None, metric='pearson'):
    if sample_list:
        sample_list = file2list(sample_list)
    if enh_state_list:
        enh_state_list = file2list(enh_state_list)
    # Load in the data pickle or process the data de-novo:
    data = get_data(data_dir, sample_dir, sample_list, saved)

    # Update sample list if using specific samples:
    if sample_list:
        data.sample_list = [os.path.basename(x) for x in sample_list]
    # list of dirs with enhancers (should be all directories)
    slist = [os.path.basename(s) for s in data.all_sample_dir_list]

    # Get the enhancers state matrices:
    # Conditioned on if there is a which_first sample to run on.
    if which_first and os.path.basename(which_first) in slist:
        idx = slist.index(os.path.basename(which_first))
        data.get_enh_state(enh_state_list, uniform=True, which_first=idx)
    else:
        print("using default state", file=sys.stderr)
        data.get_enh_state(enh_state_list, uniform=True)
    i_tss = data.gene_keys.index("tss")
    i_chrom = data.gene_keys.index("chrom")
    i_strand = data.gene_keys.index("strand")
    i_expr = data.gene_keys.index("expr")
    sys.stderr.flush()
    links = []
    bad_count = 0
    # Open file connection:
    # if type(output) == str:
        # outfile = open(output, 'a')
    # For each position bin, write out:
    for pos_bin in range(-data.flank // data.bin_size,
                         data.flank // data.bin_size + 1):
        classifier = LogisticRegression(
            penalty=PENALTY, C=REGULARIZATION, max_iter=1000)
        sys.stdout.flush()
        # train_data, test_data = create_train_test_matrix(pos_bin,data,metric)
        train_data, train_labels, train_genes, test_data, test_genes = \
            create_train_test_matrix(pos_bin, data, metric)
        print("bin %d gene loop done" % pos_bin, len(train_data),
              len(test_data), file=sys.stderr)
        sys.stdout.flush()
        sys.stderr.flush()
        if test_data is None:
            print("No test data, continuing", file=sys.stderr)
            continue
        # TODO: Had issue with fit - do we need to add ALL data?
        # Train classifier:
        classifier.fit(train_data, train_labels)
        score = classifier.score(train_data, train_labels)
        print("Score on training set:", score, file=sys.stderr)
        sys.stdout.flush()
        # Predict with classifier:
        prob_list = classifier.predict_proba(test_data)
        pos_class_idx = list(classifier.classes_).index(1)
        print("coef:", classifier.coef_)
        for i_gene, prob in zip(test_genes, prob_list):
            if prob[pos_class_idx] >= threshold / (1 + threshold):
                strand = data.gene_values[i_gene][i_strand]
                tss = data.gene_values[i_gene][i_tss] // data.bin_size
                if strand == "+":
                    item = [-pos_bin * data.bin_size, prob[pos_class_idx],
                            data.gene_values[i_gene][i_chrom],
                            data.bin_size * (tss - pos_bin),
                            data.bin_size * (tss - pos_bin + 1),
                            i_gene, data.gene_list[i_gene]]
                else:
                    item = [-pos_bin * data.bin_size, prob[pos_class_idx],
                            data.gene_values[i_gene][i_chrom],
                            data.bin_size * (tss + pos_bin),
                            data.bin_size * (tss + pos_bin + 1),
                            i_gene, data.gene_list[i_gene]]
                ratio = round(item[1] / (1 - item[1]), 4)
                string = "\t".join(map(str, [item[2], item[3], item[4],
                                             item[6], ratio, item[0]]))
                print(string + "\tITEM", file=sys.stderr)
                links.append(string)
                # for link in links:
                if type(output) == str:
                    with open(output, 'a') as outfile:
                        print("Writing to file")
                        outfile.write(string + '\n')
                else:
                    print(string, file=output)
    # outfile.close()


# Create the train/test data for a specific pos_bin (run on matrices)
def create_train_test_matrix(pos_bin, data, metric):
    print("pos_bin:", pos_bin, file=sys.stderr)
    sys.stdout.flush()
    sys.stderr.flush()
    i_tss = data.gene_keys.index("tss")
    i_chrom = data.gene_keys.index("chrom")
    i_strand = data.gene_keys.index("strand")
    i_expr = data.gene_keys.index("expr")
    # Setup classifier:
    train_data = []
    train_labels = []
    train_genes = []
    test_data = []
    test_genes = []
    ispos = []

    # For storing data:
    mark_chromlist = []
    mark_poslist = []
    expr_list = []
    rand_list = []
    rand_gene_indices = data.get_random_gene_order()
    t0 = time.time()

    for s_bin in range(pos_bin - data.smooth // data.bin_size,
                       pos_bin + data.smooth // data.bin_size + 1):
        # print("s_bin:", s_bin, pos_bin - data.smooth // data.bin_size,
        #      pos_bin + data.smooth // data.bin_size + 1, file=sys.stderr)
        sys.stdout.flush()
        for i_gene, gene in enumerate(data.gene_list):
            # print("beginning", i_gene, "in", s_bin)
            sys.stdout.flush()
            # Get expr values for the genes in order vs shuffled:
            v_gene = data.gene_values[i_gene]
            r_gene = data.gene_values[rand_gene_indices[i_gene]]
            tss_bin = v_gene[i_tss] // data.bin_size
            # all_tss[i_gene]
            # all_chrom[i_gene]
            v_chrom = v_gene[i_chrom]
            # Find position s_bin upstream of TSS (may be negative s_bin)
            if v_gene[i_strand] == "+":
                shift_pos = tss_bin - s_bin
            else:
                shift_pos = tss_bin + s_bin
            # NOTE: Skip if position is not in genome, or state is not enhancer:
            if shift_pos < 0 or \
                    shift_pos >= data.enh_state[0][v_chrom].shape[0] or \
                    data.enh_state[0][v_chrom][shift_pos] == 0:
                # must be in correct state
                continue
            # Expression values for the gene + matched
            expr_vals = v_gene[i_expr]
            # NOTE: Skip if no expression for gene:
            if np.sum(expr_vals) == 0:
                continue
            random_expr_vals = r_gene[i_expr]
            # Add to list:
            expr_list.append(expr_vals)
            rand_list.append(random_expr_vals)
            # Append indices (same across marks):
            mark_chromlist.append(v_chrom)
            mark_poslist.append(shift_pos)
            # For each gene + rand gene, add the data to training:
            train_labels.append(1)
            train_genes.append(i_gene)
            train_labels.append(0)
            train_genes.append(rand_gene_indices[i_gene])
            # Keep track of testing data pos:
            if s_bin == pos_bin:
                ispos.append(1)
            else:
                ispos.append(0)

    # Compute all of the correlations:
    t1 = time.time()
    print('Timing for index load:', str(round(t1 - t0, 2)) + 's')
    mark_chromind = {}
    for v_chrom in data.chromlist:
        mark_chromind[v_chrom] = []
    for i, val in enumerate(mark_chromlist):
        mark_chromind[val].append(i)
    # Make the datasets, calc chrom:
    mark_poslist = np.array(mark_poslist)
    expr_list = np.array(expr_list)
    rand_list = np.array(rand_list)
    pos_data = np.zeros((len(mark_poslist), len(data.marks)))
    neg_data = np.zeros((len(mark_poslist), len(data.marks)))
    train_data = np.zeros((len(mark_poslist)*2, len(data.marks)))
    pt2 = time.time()
    for j, v_mark in enumerate(data.marks):
        xfull = np.zeros((len(mark_poslist), expr_list.shape[1]))
        for chrom in data.chromlist:
            cind = np.array(mark_chromind[chrom])
            pos = mark_poslist[cind]
            x = data.mark_table_list[v_chrom][v_mark][cind]
            x = x.toarray()
            xfull[cind,:] = x

        pos_corr_vec = np.array(pearson_mat(xfull, expr_list))
        neg_corr_vec = np.array(pearson_mat(xfull, rand_list))
        pos_corr_vec[np.isnan(pos_corr_vec)] = 0
        neg_corr_vec[np.isnan(neg_corr_vec)] = 0
        pos_data[:,j] = pos_corr_vec
        neg_data[:,j] = neg_corr_vec
    pt4 = time.time()
    print(pt4 - pt2)

    # Compose the pos data and neg data into training data:
    for i in range(pos_data.shape[0]):
        train_data[(i * 2),:] = pos_data[i,:]
        train_data[(i * 2 + 1),:] = neg_data[i,:]
    # Make the test dataset where s_bin == pos_bin:
    posind = np.array(np.where(ispos)[0])
    test_data = pos_data[posind,:]
    test_genes = [train_genes[2 * pi] for pi in posind]
    t2 = time.time()
    print('Timing for corr calc:', str(round(t2 - t1, 2)) + 's')
    return(train_data, train_labels, train_genes, test_data, test_genes)


# Scans a text file into a list
def file2list(fname, split=True):
    out = []
    with open(fname, 'r') as R:
        for rline in R:
            line = rline.strip()
            if line:
                if split:
                    out += line.split()
                else:
                    out.append(line)
    return out


if __name__ == "__main__":
    fire.Fire(run)
