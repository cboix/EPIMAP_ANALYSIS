#!/usr/bin/python
# ==================================
# Turn the mtx temporary files into
# MTX files + attribute files
# Sparse matrix + attributes (pickle)
# ==================================
import time
import numpy as np
import datetime
import argparse
import csv
import gzip
import six.moves.cPickle as pickle
from scipy.sparse import csr_matrix, coo_matrix


class mtx_parser(object):
    def __init__(self, mtxfile, cellfile, outprefix):
        # rows, cols):
        # NOTE: Assumes mtx file sorted by pk
        self.mtxfile = mtxfile
        self.cellfile = cellfile
        self.outprefix = outprefix
        # Initialize the matrix information:
        self.irows = []
        self.icols = []
        self.idata = []
        self.orows = []
        self.ocols = []
        self.odata = []
        # -----------------------------------
        self.today = datetime.datetime.today().strftime("%Y%m%d")
        self.attr_file = self.outprefix + '_attr'
        # Attributes:
        self.cd = {}
        self.pd = {}

    def save_pickle(self, filename, matrix):
        out = csr_matrix(matrix)
        with open(filename, 'wb') as f:
            pickle.dump(out, f, protocol=2)

    def save_mtx(self, filename, matrix):
        coo = coo_matrix(matrix)
        with gzip.open(filename, 'wb') as f:
            # Write header:
            f.write(bytes("%%MatrixMarket matrix coordinate real general\n",
                          'UTF-8'))
            f.write(bytes(str(coo.shape[0]) + " " +
                          str(coo.shape[1]) + " " +
                          str(coo.nnz) + "\n", 'UTF-8'))
            # Write lines:
            for i in range(coo.nnz):
                f.write(bytes(str(coo.row[i]) + " " +
                              str(coo.col[i]) + " " +
                              str(coo.data[i]) + "\n", 'UTF-8'))

    # Read in cells metadata/ids:
    def load_cells(self):
        print("[STATUS] Loading column (cell) information")
        with open(self.cellfile, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for line in reader:
                self.cd[line[0]] = int(line[1])

    # Read in matrix line by line:
    def process_mat(self):
        self.start = time.time()
        self.rid = -1
        self.rid_last = -1
        self.pk = ""
        print("[STATUS] Reading in lines")
        with open(self.mtxfile, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for line in reader:
                self.process_line(line)
        print("[STATUS] Making matrices")
        self.make_matrices()

    # Process + save lines to COO matrices:
    def process_line(self, line):
        if line[0] != self.pk:
            self.pk = line[0]
            if self.pk in self.pd.keys():
                self.rid = self.pd[self.pk]
            else:
                self.rid = self.rid_last + 1
                self.rid_last = self.rid
                self.pd[self.pk] = self.rid
            # Print progress.
            if self.rid % 100000 == 0:
                print("[STATUS] Line # " + str(self.rid) +
                      " TIME: " + str(round(time.time() - self.start, 2)) + "s")
        self.cid = self.cd[line[2]] - 1  # Change to 0 indexing.
        # NOTE: Ignores possibility of multi-pk from same ct, but we throw these
        # large regions out afterwards.
        if line[1] == 'imputed':
            self.irows.append(self.rid)
            self.icols.append(self.cid)
            self.idata.append(float(line[3]))
        else:
            self.orows.append(self.rid)
            self.ocols.append(self.cid)
            self.odata.append(float(line[3]))

    # Make matrices from IJV formatted lines:
    def make_matrices(self):
        M = len(self.pd)
        N = len(self.cd)
        print(M)
        print(self.rid_last)
        self.imp = coo_matrix((np.array(self.idata).astype(np.int),
                               (np.array(self.irows).astype(np.int),
                                np.array(self.icols).astype(np.float))),
                              shape=(M, N))
        self.obs = coo_matrix((np.array(self.odata).astype(np.int),
                               (np.array(self.orows).astype(np.int),
                                np.array(self.ocols).astype(np.float))),
                              shape=(M, N))

    def save_matrices(self):
        # Save sparse matrices:
        print("[STATUS] Saving imputed CSR")
        self.save_pickle(self.outprefix + "_imputed_csr.cp", self.imp)
        print("[STATUS] Saving observed CSR")
        self.save_pickle(self.outprefix + "_observed_csr.cp", self.obs)
        # Save mtx data:
        print("[STATUS] Saving imputed MTX")
        self.save_mtx(self.outprefix + "_imputed.mtx.gz", self.imp)
        print("[STATUS] Saving observed MTX")
        self.save_mtx(self.outprefix + "_observed.mtx.gz", self.obs)

    # Save the dictionaries:
    def save_attributes(self):
        print("[STATUS] Saving the attributes (peaks + cells)")
        self.save_attr_cp()
        # Peaks to TSV:
        with gzip.open(self.outprefix + "_rows.tsv.gz", 'wb') as f:
            for key in self.pd:
                f.write(bytes(str(self.pd[key]) + "\t" + str(key) + "\n",
                              'UTF-8'))
        # Cells to TSV:
        with gzip.open(self.outprefix + "_cols.tsv.gz", 'wb') as f:
            for key in self.cd:
                f.write(bytes(str(self.cd[key]) + "\t" + str(key) + "\n",
                              'UTF-8'))

    def save_attr_cp(self):
        attr = {'peaks': self.pd,
                'cells': self.cd}
        with open(self.outprefix + "_attr.cp", 'wb') as f:
            pickle.dump(attr, f, protocol=2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Turn txt bed files into one sparse file")
    parser.add_argument('--mtx', metavar='mtx', type=str,
                        help='temp mtx file')
    parser.add_argument('--cells', metavar='cells', type=str,
                        help='List of cells + cellid')
    parser.add_argument('--out', metavar='out', type=str,
                        help='Output prefix (with directory)')
    # Size of output matrix:
    # parser.add_argument('--nrow', metavar='nrow', type=str,
    #                     help='Number of rows (locations)')
    # parser.add_argument('--ncol', metavar='ncol', type=str,
    #                     help='Number of columns (cells)')
    # Training parameters:
    args = parser.parse_args()
    parser = mtx_parser(mtxfile=args.mtx,
                        cellfile=args.cells,
                        outprefix=args.out)
    # rows=args.nrow,
    # cols=args.ncol)
    print("---------------------------------------------------")
    parser.load_cells()
    print("---------------------------------------------------")
    parser.process_mat()
    print("---------------------------------------------------")
    parser.save_attributes()
    print("---------------------------------------------------")
    parser.save_matrices()
