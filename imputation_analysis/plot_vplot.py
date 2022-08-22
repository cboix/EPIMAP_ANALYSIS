#!/usr/bin/env python
# Author: Jason Buenrostro, Stanford University
# modified from plotV_vC.py (Alicia)
# Will make a V-plot from bed regions
import os
import sys
import numpy as np
import pysam
import matplotlib
from multiprocessing import Pool
from optparse import OptionParser
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# ------------------------------
# Read options from command line
# ------------------------------
opts = OptionParser()
usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-a", help="<Reads> Accepts sorted BAM file")
opts.add_option("-b", help="<Bed> Accepts bed file")
opts.add_option("-o", help="OutputFile")
opts.add_option("-e", default="1000",
                help="number of bases to extend to each side, default=1000")
opts.add_option("-p", default="center",
                help="options:center,ends, default=center")
opts.add_option("-c", default="20", help="number of threads to use, default=20")
opts.add_option("-r", default="", help="replace either 'chr' or '' in bedfile")
opts.add_option("-l", default="-99", help="set fixed template length")
opts.add_option(
    "-s", default='4',
    help="column in which strand information is located (1=first), default=4")
opts.add_option("-u", action="store_true", default=False,
                help="Print uncompressed output file")
opts.add_option("-v", action="store_true", default=False,
                help="Print profile around bed file")
opts.add_option("-i", action="store_true", default=False,
                help="Print insert sizes across intervals")
opts.add_option("-x", action="store_true", default=False,
                help="Plot and print interval examples")
opts.add_option("--window", default='20', help="window size for plotting")
options, arguments = opts.parse_args()

# Return usage information if no argvs given
if len(sys.argv) == 1:
    os.system(sys.argv[0] + " --help")
    sys.exit()

read_length = int(options.l)


# ----------
# Functions:
# ----------
# Assign matrix (only marks the center position!)
def asn_mat(center_position, mat, start_interval, end_interval, t, i, weight):
    if float(center_position) >= start_interval and \
            float(center_position) < end_interval - 1 and \
            t < num_isize:
        base = center_position - start_interval
        if len(intervals[0]) == 3:
            mat[t, base] += weight
        elif intervals[i][int(options.s) - 1] == "-":
            mat[t, len(mat[0]) - base - 1] += weight
        else:
            mat[t, base] += weight
    return mat


def add_read(mat, left_position, start_interval, end_interval,
             insert_length, i, full=False):
    if full:
        t = i
    else:
        t = insert_length
    right_position = left_position + insert_length
    center_position = left_position + insert_length / 2
    if insert_length % 2 == 1 and options.p == 'center':
        # Even split of read:
        mat = asn_mat(center_position, mat,
                      start_interval, end_interval, t, i, 0.5)
        mat = asn_mat(center_position + 1, mat,
                      start_interval, end_interval, t, i, 0.5)
    elif insert_length % 2 != 1 and options.p == 'center':
        # By insert-size:
        mat = asn_mat(center_position, mat,
                      start_interval, end_interval, t, i, 1)
    elif options.p == 'ends':
        # Add read to either end.
        mat = asn_mat(left_position, mat,
                      start_interval, end_interval, t, i, 1)
        mat = asn_mat(right_position, mat,
                      start_interval, end_interval, t, i, 1)
    else:
        sys.exit('Error, check parameters')
    return(mat)


# Compute vPlot Matrix for a particular chunk:
def sub_Mat(start):
    # Initialize data matrix
    mat = np.zeros([num_isize, num_bp])
    intmat = np.zeros([num_intervals, num_bp])
    # Loop through the intervals and get relevant info
    bamfile = pysam.Samfile(options.a, "rb")
    end = min(start + chunksize, len(intervals))
    for i in range(start, end):
        # Get center, start, and end of interval (extend from center)
        center = int(intervals[i][1]) + \
            (int(intervals[i][2]) - int(intervals[i][1])) / 2
        start_interval = center - int(options.e)
        end_interval = center + int(options.e)
        # Loop through reads
        for reads in bamfile.fetch(str(intervals[i][0]),
                                   max(0, start_interval - 2000),
                                   end_interval + 2000):
            # Filter on mapping quality
            if reads.mapq < 30:  # or reads.is_proper_pair==False:
                continue
            # Filter on reversed:
            if reads.is_reverse:
                continue
            else:
                # Calculate read positions:
                readshift = 0  # Trimming `readshift` bp off each side
                left_position = reads.pos + readshift
                # Calculate center point
                if read_length > 0:
                    insert_length = abs(read_length) - 1 - readshift * 2
                else:
                    insert_length = abs(reads.template_length) - 1 - \
                        readshift * 2
                # Add to insert-size matrix:
                mat = add_read(mat, left_position, start_interval,
                               end_interval, insert_length, i, full=False)
                # Add to intervals matrix:
                intmat = add_read(intmat, left_position, start_interval,
                                  end_interval, insert_length, i, full=True)
    return(mat, intmat)

# -----------
# Run script:
# -----------
# Get intervals and change contig to non-chr labeled
intervals = np.loadtxt(options.b, 'str')
for i in range(len(intervals)):
    intervals[i][0] = intervals[i][0].replace(options.r, "")

# Determine number of rows and columns for matrix
num_isize = 1000  # Number of insert sizes.
num_bp = int(options.e) * 2  # Width of region in bps (direct 2x)
# Extension:
# num_bp = int(intervals[0][2]) - int(intervals[0][1]) + int(options.e)*2

# Split bedfile into chunks
num_intervals = len(intervals)
chunksize = num_intervals / int(options.c)
chunks = num_intervals / chunksize
starts = range(0, num_intervals, chunksize)
print("Chunked into: " + str(starts))

# Computation of matrix:
if int(options.c) == 1:
    # Single-core:
    (submat, subint) = sub_Mat(0)
    mat = submat
    intmat = subint
else:
    # Parallel processing for each chunk:
    if __name__ == "__main__":
        pool = Pool(processes=int(options.c))
        # TODO fix this for int/sub mat for multiple.
        sub_mats = pool.map(sub_Mat, starts, 1)
    # Sum up matrices for each chunk into matrix for all:
    mat = np.zeros([num_isize, num_bp])
    for i in range(len(starts)):
        mat = mat + sub_mats[i]

# Get margins (vplot and insert-size):
vmat = np.sum(mat, 0)
imat = np.sum(mat, 1)

# =============================
# Save matrices and make plots:
# =============================
# Prefixes for files:
if not options.o:
    n1 = os.path.basename(options.a)
    n2 = os.path.basename(options.b)
    basename = n1 + "." + n2
else:
    basename = options.o
print("Using basename: " + basename)

# Save matrix:
ufile = basename + '.vplot'
vfile = basename + '.vect'
ifile = basename + '.iSize'
xfile = basename + '.examples'
if options.u:
    np.savetxt(ufile + '.tsv', mat, delimiter='\t', fmt='%s')
if options.v:
    np.save(vfile, vmat)
if options.i:
    np.save(ifile, imat)
if options.x:
    np.savez_compressed(xfile, intmat)

# Plot vplot profile
xran = min(500, int(options.e))
yran = min(500, num_isize)
if options.v:
    print("Plotting v-plot")
    fig = plt.figure(figsize=(8.0, 5.0))
    ax = fig.gca()
    # plt.plot(mat/np.median(mat[1:200]))
    plt.plot(vmat / np.mean(vmat[0:yran]), 'k.')
    plt.plot(np.convolve(vmat, np.ones(int(options.window)), 'same') /
             int(options.window) / np.mean(vmat[0:yran]), 'r')
    plt.xlabel('Position relative to center')
    plt.ylabel('Insertions')
    axvline(xran, color='grey', ls='dashed')
    fig.savefig(vfile + '.png')

# Plot insert-sizes:
if options.i:
    print("Plotting insert-sizes")
    fig = plt.figure(figsize=(8.0, 5.0))
    plt.plot(imat[0:990])
    plt.xlabel('Insert size (bp)')
    plt.ylabel('Count')
    fig.savefig(ifile + '.png')

# Plot heatmap of profile:
if options.u:
    print("Plotting heatmap of vplot at different insert sizes")
    fig = plt.figure(figsize=(8.0, 5.0))
    xstart = int(options.e) - xran
    xend = int(options.e) + xran + 1
    print(0, yran)
    print(xstart, xend)
    plt.imshow(mat[0:yran, xstart:xend],
               origin='lower', aspect='equal',
               extent=[-xran, xran, 1, yran + 1])
    plt.xlabel('Position relative to center')
    plt.ylabel('Insert size')
    fig.savefig(ufile + '.png')

# Sort intmat by total number of counts:
print(intmat.shape)
intcounts = np.sum(intmat, 1)
print(len(intcounts))
order = np.argsort(-intcounts, 0)
print(len(order))
intsorted = intmat[order, :]
nexamples = 2000

# Plot heatmap of examples:
if options.u:
    print("Plotting heatmap of examples")
    fig = plt.figure(figsize=(7.0, 11.0))
    ax = fig.gca()
    xstart = int(options.e) - xran
    xend = int(options.e) + xran + 1
    pltmat = np.log2(intsorted[0:nexamples, xstart:xend] + 1)
    plt.imshow(pltmat, vmin=0, vmax=3, cmap='Blues',
               extent=[-xran, xran, nexamples + 1, 0])
    ax.axvline(0, color='grey', ls='dashed')
    plt.xlabel('Position relative to center')
    plt.ylabel('Interval number')
    fig.savefig(xfile + '.png')
