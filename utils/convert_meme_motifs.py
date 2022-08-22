#!/usr/bin/pyton
# ------------------------------
# Utility to convert MEME motifs
# to pfm or to Pouya's format
# ------------------------------
# cmd="python $BINDIR/convert_meme_motifs.py main --motiffile jolma2013.meme --outprefix jolma2013 --nomouse"
# cmd="python $BINDIR/convert_meme_motifs.py main --motiffile JASPAR2018_CORE_vertebrates_non-redundant_pfms.meme --outprefix JASPAR2018_CORE_vertebrates_non-redundant_pfms --nomouse"
# cmd="python $BINDIR/convert_meme_motifs.py main --motiffile HOCOMOCOv11_core_HUMAN_mono_meme_format.meme --outprefix HOCOMOCOv11"
import re
import os
import fire
import numpy as np


class convert_meme_motifs(object):
    def __init__(self, motiffile, outprefix, nomouse=False):
        self.motiffile = motiffile
        self.outprefix = outprefix
        self.nomouse = nomouse

    def main(self):
        self.read_motifs()
        self.write_pfm()
        self.write_motifs()

    def read_motifs(self):
        print("[STATUS] Reading motifs")
        self.names = []
        self.attr = []
        self.pfm = []
        inmotif = 0
        w = 0
        with open(self.motiffile, 'r') as f:
            for line in f:
                line = line.rstrip("\n")
                line = line.split(" ")
                if line[0] == "MOTIF":
                    self.names.append(line[1])
                    mat = []
                    print(line[1])
                elif w > 0:
                    w = w - 1
                    mat.append([v for v in line if v != ''])
                    # mat.append([float(v) for v in line if v != ''])
                    if w == 0:
                        self.pfm.append(mat) # Append motif
                elif line[0] == "letter-probability":
                    inmotif = 1
                    self.attr.append(" ".join(line))
                    # Extract width of motif:
                    wind = [i+1 for i,v in enumerate(line) if v == 'w='][0]
                    w = int(line[wind])
                    print('Length is', w)
                elif line != ['']:
                    print(line)

    def write_pfm(self):
        print("[STATUS] Writing PFM format")
        with open(self.outprefix + "_pfm.txt", 'w') as f:
            if self.nomouse:
                mlist = [i for i, v in enumerate(self.names) 
                        if not re.search('mouse', v)]
            else:
               mlist = range(len(self.names))
            for i in mlist:
                f.write(">" + str(self.names[i]) + " " + self.attr[i] +  "\n")
                mat = np.array(self.pfm[i]).T
                mat = mat.tolist()
                for j in range(len(mat)):
                    f.write(" ".join(mat[j]) + "\n")

    def write_motifs(self):
        print("[STATUS] Writing motifs.txt format")
        with open(self.outprefix + "_motifs.txt", 'w') as f:
            if self.nomouse:
                mlist = [i for i, v in enumerate(self.names) 
                        if not re.search('mouse', v)]
            else:
               mlist = range(len(self.names))
            for i in mlist:
                f.write(">" + str(self.names[i]) + " " + self.attr[i] +  "\n")
                for j in range(len(self.pfm[i])):
                    f.write("X " + " ".join(self.pfm[i][j]) + "\n")

if __name__ == "__main__":
    fire.Fire(convert_meme_motifs)
