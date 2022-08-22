#!/usr/bin/python
# ----------------------------------------------
# Utility to calculate difference between tracks
# Rescales if prompted
# Outputs slope and plots
# ----------------------------------------------
import os
import re
import gzip
import socket
import time
import numpy as np

# For cli usage:
import fire

rng = np.random


class diff_bws(object):
    def __init__(self, observed, imputed, output, rescale=True):
        # Arguments:
        self.observed = observed
        self.imputed = imputed
        self.output = output
        self.rescale = rescale
        # Check domain and set home:
        self.domain = socket.getfqdn()
        if 'broadinstitute.org' in self.domain:
            self.HOME = '/broad/compbio/cboix'
        else:
            self.HOME = os.getenv("HOME")
        self.datadir = self.HOME + '/EPIMAP_ANALYSIS/db/'
        self.trackdir = self.datadir + 'ChromImpute/imputed/'
        self.statsdir = self.datadir + 'ChromImpute/stats/'
        self.chrlist = ['chr' + str(i+1) for i in range(22)] + ['chrX']
        self.start = time.time()
        # Initialize regression parameters:
        self.X_avg = 0
        self.Y_avg = 0
        self.Sxy = 0
        self.Sx = 0
        self.n = 0
        self.alpha = 0
        self.beta = 0

    def main(self):
        # Perform regression online:
        self.perform_regression()
        print("Regression took: " +
              str(round(time.time() - self.start, 2)) + "s")
        # Write output of calculated slope:
        self.write_slopes()
        # Re-read chromosomes and output gzipped difference files:
        self.write_differences()
        print("Regression + writing took: " +
              str(round(time.time() - self.start, 2)) + "s")

    def perform_regression(self):
        for chrom in self.chrlist:
            print(str(chrom) + ". Slope is " + str(round(self.beta, 2)))
            # - Read in imp/obs
            prefix = self.trackdir + chrom + "_"
            [self.obs, h1, h2] = read_bw(prefix + self.observed + ".wig.gz")
            [self.imp, h1, h2] = read_bw(prefix + self.imputed + ".wig.gz")
            # - Update regression coefficients
            [self.Sxy, self.Sx,
             self.n, self.alpha, self.beta,
             self.X_avg, self.Y_avg] = update_lr(self.X_avg, self.Y_avg,
                                                 self.Sxy, self.Sx, self.n,
                                                 self.obs, self.imp)

    def write_differences(self):
        for chrom in self.chrlist:
            print("Re-loading and writing: " + str(chrom))
            # - Read in imp/obs
            prefix = self.trackdir + chrom + "_"
            [self.obs, h1, h2] = read_bw(prefix + self.observed + ".wig.gz")
            [self.imp, h1, h2] = read_bw(prefix + self.imputed + ".wig.gz")
            # - Calculate obs - imp / slope, and clamp values at 0 (for corr.)
            self.diff = np.clip(self.obs - self.imp / self.beta,
                    a_min=0, a_max=None)
            # - Change the header lines to reflect diff
            h1diff = re.sub("_imputed", "_difference", h1)
            # - Write out the difference using imputed prefix
            write_bw(self.diff, [h1diff, h2], prefix + self.output + ".wig.gz")

    def write_slopes(self):
        # Write coefficients to the stats directory:
        with open(self.statsdir + self.output + "_reg.tsv", "w") as f:
            f.write(str(round(self.alpha, 4)) + "\t" +
                    str(round(self.beta, 4)) + "\n")


# Online linear regression, as in:
# https://stackoverflow.com/questions/52070293/
def update_lr(x_avg, y_avg, Sxy, Sx, n, new_x, new_y):
    """
    x_avg: average of previous x, if no previous sample, set to 0
    y_avg: average of previous y, if no previous sample, set to 0
    Sxy: covariance of previous x and y, if no previous sample, set to 0
    Sx: variance of previous x, if no previous sample, set to 0
    n: number of previous samples
    new_x: new incoming 1-D numpy array x
    new_y: new incoming 1-D numpy array x
    """
    # Update stats:
    new_n = n + len(new_x)
    new_x_avg = (x_avg * n + np.sum(new_x)) / new_n
    new_y_avg = (y_avg * n + np.sum(new_y)) / new_n
    if n > 0:
        x_star = (x_avg * np.sqrt(n) + new_x_avg * np.sqrt(new_n)) / \
            (np.sqrt(n) + np.sqrt(new_n))
        y_star = (y_avg * np.sqrt(n) + new_y_avg * np.sqrt(new_n)) / \
            (np.sqrt(n) + np.sqrt(new_n))
    elif n == 0:
        x_star = new_x_avg
        y_star = new_y_avg
    else:
        raise ValueError
    # Update variance and covariance:
    new_Sx = Sx + np.sum((new_x - x_star)**2)
    new_Sxy = Sxy + np.sum((new_x - x_star).reshape(-1) *
                           (new_y - y_star).reshape(-1))
    # Update coefficients
    beta = new_Sxy / new_Sx
    alpha = new_y_avg - beta * new_x_avg
    return(new_Sxy, new_Sx, new_n, alpha, beta, new_x_avg, new_y_avg)


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


def write_bw(array, header, file):
    with gzip.open(file, 'wb') as f:
        # Write two header lines:
        f.write(header[0])
        f.write(header[1])
        for item in array:
            f.write(str(round(item, 3)) + "\n")

if __name__ == "__main__":
    fire.Fire(diff_bws)
