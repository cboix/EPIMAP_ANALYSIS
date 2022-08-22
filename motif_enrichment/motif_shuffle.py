#!/usr/bin/python
# --------------------------------------------
# Make shuffled motifs from a pwm or pfm file:
# --------------------------------------------
import os
import re
import fire
import numpy as np

rng = np.random


class shuffle_motif(object):
    def __init__(self, infile, outprefix, nshuffle=25, seed=1):
        # Arguments:
        self.infile = infile
        self.outprefix = outprefix
        self.nshuffle = nshuffle
        self.seed = seed
        self.name = None
        rng.seed(self.seed)

    def main(self):
        self.read_motif()
        self.write_shuffles()

    def read_motif(self):
        print("Reading motif")
        self.mat = []
        inmat = False
        with open(self.infile, 'r') as f:
            for line in f:
                line = line.rstrip("\n")
                if re.search("^>", line):
                    self.name = re.sub("^>","",line).split(" ")
                elif line != [""]:
                    line = line.split()
                    line = [float(v) for v in line]
                    self.mat.append(line)
        # Determine axis for shuffling:
        self.mat = np.array(self.mat)

    def write_shuffles(self):
        print("Writing shuffles")
        for i in range(self.nshuffle):
            # Shuffle:
            smat = self.mat
            if self.mat.shape[0] == 4:
                w = self.mat.shape[1]
                ind = np.array(rng.choice(range(w), size=w, replace=False))
                smat = smat[:,ind]
            elif self.mat.shape[1] == 4:
                w = self.mat.shape[0]
                ind = np.array(rng.choice(range(w), size=w, replace=False))
                smat = smat[ind,:]
            else:
                print("Axis was neither 0 nor 1.")
                break
            # Write the shuffles out:
            fname = self.outprefix + "_shuffled_" + str(i) + '.txt'
            with open(fname, 'w') as f:
                smat = smat.tolist()
                if self.name is not None:
                    self.name[0] = self.name[0] + "_shuffled_" + str(i)
                    f.write(">" + " ".join(self.name) + "\n")
                for line in smat:
                    line = [str(v) for v in line]
                    line = " ".join(line) + "\n"
                    f.write(line)


if __name__ == "__main__":
    fire.Fire(shuffle_motif)
