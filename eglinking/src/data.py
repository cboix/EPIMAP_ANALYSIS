import numpy as np
import os
import sys
from scipy.sparse import csr_matrix, coo_matrix, csc_matrix
from bed import get_bed, get_chrom_size_table

CHROM_FILE = "hg19.chrom.sizes"
EXPR_FILE = "rna.txt"
TSS_FILE = "tss.txt"

class Data:
    def __init__(self, data_dir, sample_outer_dir, sample_list, mark_list,
                 bin_size=200, flank=1000000, smooth=5000, mask=None):
        # 0. Parameters + arguments
        self.data_dir = data_dir
        self.sample_dir_list = [os.path.join(sample_outer_dir, s)
                                for s in sample_list]
        # in case dirnames were also provided
        self.sample_list = [os.path.basename(s)
                            for s in sample_list]
        self.all_sample_dir_list = list(
            filter(os.path.isdir, [os.path.join(sample_outer_dir, s)
                                   for s in os.listdir(sample_outer_dir)
                                   if s[:3] == 'BSS']))
        self.nsamp = len(self.sample_dir_list)
        self.mask_regions = None
        self.bin_size = bin_size
        self.flank = flank
        self.smooth = smooth
        self.chromlist = None
        self.permark = True  # Create a per-mark dict
        # 1. Get gene list for the dataset:
        gene_dir_list = [os.path.join(sample_dir, EXPR_FILE)
                         for sample_dir in self.sample_dir_list]
        self.get_gene_list(gene_dir_list,
                     os.path.join(data_dir, TSS_FILE))
        chrom_file = os.path.join(self.data_dir, CHROM_FILE)
        # 2. Load genome information + mask if used
        if not os.path.isfile(chrom_file):
            print("%s is not a valid file" % chrom_file, file=sys.stderr)
            sys.exit(1)
        self.chrom_size = get_chrom_size_table(chrom_file)
        # Training mask (e.g. if we only want to work with chr1):
        if mask:
            with open(mask, "r") as R:
                self.mask_regions = get_bed(self.chrom_size, self.bin_size, R)
        # 3. Load data across samples for all marks:
        self.get_mark(mark_list)

    def get_mark(self, mark_list):
        self.marks = sorted(list(mark_list))
        self.mark_table_list = []
        mtl = {}
        for i, sample_dir in enumerate(self.sample_dir_list):
            mark_table = {}
            print("reading marks for", os.path.basename(sample_dir),
                  file=sys.stderr)
            sys.stdout.flush()
            fileok = {}
            for j, mark in enumerate(self.marks):
                print(j, "- mark:", mark)
                mark_file = os.path.join(sample_dir, mark)
                if not os.path.isfile(mark_file):
                    print("%s is not a valid file" % mark_file, file=sys.stderr)
                    fileok[mark] = False
                    # sys.exit(1)
                else:
                    fileok[mark] = True
                    print("Reading bedfile: " + str(mark_file))
                    with open(mark_file, "r") as R:
                        # TODO: Could replace w/ sparse mat for load and write.
                        mark_bed = get_bed(self.chrom_size, self.bin_size, R)
                if self.chromlist is None:
                    self.chromlist = list(mark_bed.keys())
                    print(self.chromlist)
                if fileok[mark]:
                    if self.permark:
                        print("Concatenating with other marks")
                        print("Adding", mark, "to dictionary")
                        for chrom in self.chromlist:
                            if j == 0:
                                mark_table[chrom] = {}
                            mark_table[chrom][mark] = csc_matrix(
                                mark_bed[chrom]).transpose()
                    else:
                        if not mark_table:
                            mark_table = mark_bed
                        else:
                            print("Concatenating with other marks")
                            for chrom in mark_table.keys():
                                # concatenate arrays in mark direction
                                mark_table[chrom] = np.column_stack(
                                    (mark_table[chrom], mark_bed[chrom]))
                        # mtl[chrom][position][index]
            if i == 0:
                print("Allocating memory")
                print(mark_table.keys())
                for chrom in mark_table.keys():
                    if self.mask_regions and sum(self.mask_regions[chrom]) == 0:
                        continue
                    if self.permark:
                        for j, mark in enumerate(self.marks):
                            print(chrom, mark)
                            shp = mark_table[chrom][mark].shape
                            print(shp)
                            # shp = mark_table[chrom].shape
                            print("allocating ", shp[0], 'by',
                                  self.nsamp, file=sys.stderr)
                            sys.stdout.flush()
                            if j == 0:
                                mtl[chrom] = {}
                            # TODO: Change to lil_matrix?
                            mtl[chrom][mark] = csc_matrix((shp[0], self.nsamp))
                    else:
                        shp = mark_table[chrom].shape
                        print("allocating ", shp, 'by',
                            len(self.sample_dir_list), file=sys.stderr)
                        sys.stdout.flush()
                        mtl[chrom] = np.zeros((shp[0], shp[1], self.nsamp))
            # mtl indices: chrom -> bin -> mark -> sample
            print("Adding to mtl")
            for chrom in mtl.keys():
                # Add each mark + chrom combination at sample i
                if self.permark:
                    for j, mark in enumerate(self.marks):
                        if fileok[mark]:
                            mtl[chrom][mark][:, i] = mark_table[chrom][mark]
                else:
                    print(mark_table[chrom].shape)
                    if fileok[mark]:
                        mtl[chrom][:, :, i] = mark_table[chrom]
        print("Finished reading in marks", file=sys.stderr)
        sys.stdout.flush()
        self.mark_table_list = mtl


    def get_enh_state(self, enh_state_list, uniform=True, which_first=None):
        self.enh_state = []
        enh_bed_list_list = [[os.path.join(sample_dir, "%s.bed" % s)
                              for s in enh_state_list]
                             for sample_dir in self.all_sample_dir_list]
        if which_first and which_first in range(len(enh_bed_list_list)):
            enh_bed_list_list = [enh_bed_list_list[which_first]] + \
                enh_bed_list_list[1:]
        else:
            which_first = self.all_sample_dir_list[0]
        print("Using ChromHMM states (", ", ".join(enh_state_list),
              ") for sample ", which_first, file=sys.stderr)
        if not uniform:
            print("caution: implemented for uniform-state only",
                  file=sys.stderr)
            sys.stderr.flush()
        for enh_bed_list in enh_bed_list_list:
            total = None
            if uniform and len(self.enh_state) > 0:
                self.enh_state.append(self.enh_state[0])
                continue
            for f in enh_bed_list:
                if not os.path.isfile(f):
                    print("%s is not a valid file" % f, file=sys.stderr)
                    sys.exit(1)
                with open(f, 'r') as R:
                    res = get_bed(self.chrom_size, self.bin_size, R)
                    if not total:
                        total = res
                    else:
                        for chrom, val in res.items():
                            total[chrom] += val
                if self.mask_regions: # if we are masking any locations
                    for chrom in self.mask_regions.keys():
                        total[chrom] *= self.mask_regions[chrom]
                for chrom in total.keys():
#                    print("total[%s].shape" % chrom, total[chrom].shape)
                    total[chrom][total[chrom] != 0] = 1
            self.enh_state.append(total)

    def get_gene_list(self, expr_txt_list, tss_txt):
        print("[STATUS] Loading gene list")
        if not os.path.isfile(tss_txt):
            print("%s is not a valid file" % tss_txt, file=sys.stderr)
            sys.exit(1)
        gene_index = set() # get list of names across all samples
        gene_expr = {}     # hash of samples -> names -> expr
        for expr_txt in expr_txt_list:
            # print("[STATUS] Reading thru:", expr_txt)
            if not os.path.isfile(expr_txt):
                print("%s is not a valid file" % expr_txt, file=sys.stderr)
                sys.exit(1)
            gene_expr[expr_txt] = {}
            with open(expr_txt, 'r') as ET:
                for rline in ET:
                    L = rline.strip().split()
                    name = L[0]
                    expr = float(L[1])
                    if name not in gene_index:
                        gene_index.add(name)
                    gene_expr[expr_txt][name] = expr
        print("[STATUS] Read in", len(gene_index), "genes.")
        self.gene_keys = ["chrom", "tss", "strand", "expr"]
        # Make sure no two genes have the same location, strand, expr
        unique_gene_loc = set()
        self.gene_values = []
        self.gene_list = []
        with open(tss_txt, 'r') as TT:
            for rline in TT:
                L = rline.strip().split()
                name = L[0]
                chrom = L[1]
                tss = int(L[2])
                strand = L[3]
                if name not in gene_index: # Gene may not be measured in tissue
                    continue
                expr = []
                for k in expr_txt_list:
                    if name in gene_expr[k]:
                        expr.append(gene_expr[k][name])
                    else:
                        expr.append(0)
                key = tuple(expr + [chrom, tss // self.bin_size, strand])
                if key not in unique_gene_loc:
                    value = (chrom, tss, strand, np.asarray(expr))
                    self.gene_values.append(value)
                    self.gene_list.append(name)
                    unique_gene_loc.add(key)
        print("[STATUS] Reduced to", len(self.gene_list), "genes.")

    def get_random_gene_order(self):
        indices = list(range(len(self.gene_list)))
        # which index in tuple is chrom
        tup_chrom_idx = self.gene_keys.index("chrom")
        i_chrom_list = [self.gene_values[i][tup_chrom_idx] for i in indices]
        good = {chrom: set() for chrom in set(i_chrom_list)}
        # TODO: Should we set a seed here?
        np.random.shuffle(indices)
        for i in range(len(indices)):
            orig_chrom = i_chrom_list[i]
            shuf_chrom = i_chrom_list[indices[i]]
            if orig_chrom == shuf_chrom and len(good[orig_chrom]) == 0:
                good_indices = [i for i, chrom in enumerate(i_chrom_list[i:])
                                if chrom != orig_chrom]
                new_idx = np.random.choice(good_indices)
                temp = indices[i]
                indices[i] = indices[new_idx]
                indices[new_idx] = temp
            elif orig_chrom == shuf_chrom:
                val = np.random.choice(tuple(good[orig_chrom]))
                for chrom in good.keys():
                    if val in good[chrom]:
                        good[chrom].remove(val)
                temp = indices[i]
                indices[i] = indices[val]
                indices[val] = temp
            if orig_chrom != shuf_chrom:
                for chrom in good.keys():
                    if chrom not in (i_chrom_list[indices[i]], orig_chrom):
                        good[chrom].add(i) # indices[i]
        return indices
