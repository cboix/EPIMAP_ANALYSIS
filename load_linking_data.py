#!/usr/bin/python
# --------------------------------
# Load datasets for learning links
# --------------------------------
import os
import re
import fire
import glob
import time
import gzip
import six.moves.cPickle as pickle
from scipy.sparse import coo_matrix

# For genomic region overlaps:
import parse_config
import auxfunc_masterlist as AUX
import pyranges as pr
import ncls
from pyranges.multithreaded import pyrange_apply

# For plotting points
import pandas as pd
import numpy as np

# For random data:
rng = np.random


# Load X from file or give obj
class load_linking_data(object):
    def __init__(self):
        # Directories:
        self.datadir = os.environ['HOME'] + os.environ['DBDIR']
        self.anndir = self.datadir + "/Annotation"
        self.qtldir = self.datadir + '/GTEx_eQTL/GTEx_Analysis_v7_eQTL'
        self.keptlist = self.anndir + '/kept_bssid_20190322.txt'
        self.orderfile = self.anndir + '/bssid_order_frz20190326_samponly.tsv'
        self.dmlpref = os.environ['HOME'] + os.environ['DML_DIR'] + \
            '/masterlist_DHSs_733samples_WM20180608_all_coords'
        self.re_mid = '_r200_e0'
        self.coordfile = self.dmlpref + '_hg19.core.srt.txt'
        self.coordnamfile = self.dmlpref + '_hg19' + \
            self.re_mid + "_names.core.srt.tsv"
        self.start = time.time()
        self.nclust = 300

    def main(self):
        self.load_data()

    def load_data(self):
        self.load_matrices()
        self.load_annotations()
        self.load_GTEx_eQTLs()

    # Load data matrices and element indices:
    def load_matrices(self):
        print("[STATUS] Loading matrices x5")
        # Establish dataset prefixes:
        model = 'observed_aux_18_on_mixed_impobs_QCUT_'
        mtail = '_bin_on_mixed_impobs_r200_e0_allchr_merged'
        htail = '_all_bin_on_mixed_impobs_r25_e100_allchr_merged'
        self.ecpref = self.datadir + '/' + model + 'ENH' + mtail
        self.pcpref = self.datadir + '/' + model + 'PROM' + mtail
        h1pref = self.datadir + '/' + 'H3K27ac' + htail
        h2pref = self.datadir + '/' + 'H3K4me1' + htail
        h3pref = self.datadir + '/' + 'H3K4me3' + htail
        # Load data/names:
        self.eX = get_data(self.ecpref)
        self.pX = get_data(self.pcpref)
        self.full_emarg = np.array(np.sum(self.eX, axis=1).T)[0]
        self.full_pmarg = np.array(np.sum(self.pX, axis=1).T)[0]
        self.enam = get_names(self.ecpref)
        self.pnam = get_names(self.pcpref)
        # Load + reorder columns of histone matrices:
        self.h1X = self.load_histonemat(h1pref, self.enam, binary=False)
        self.h2X = self.load_histonemat(h2pref, self.enam, binary=False)
        self.h3X = self.load_histonemat(h3pref, self.enam, binary=False)
        # Get enh/promoter indices:
        eindpref = self.datadir + '/ENH_masterlist_indices'
        pindpref = self.datadir + '/PROM_masterlist_indices'
        self.enhind = AUX.load_pickle_gzip(eindpref + "_ind.cp.gz")
        self.promind = AUX.load_pickle_gzip(pindpref + "_ind.cp.gz")

    # Load histone matrices to conform with others:
    def load_histonemat(self, prefix, nam, binary=False):
        # Load and reorder names:
        hnam = get_names(prefix)
        hnam = [re.sub(".*BSS", "BSS", n).split("_")[0] for n in hnam]
        hord = np.array([np.where(np.array(hnam) == n)[0][0] for n in nam])
        # Load and reorder matrix:
        tmpX = get_data(prefix)
        if binary:
            hXord = 1 * (tmpX[:, hord] > 0)
        else:
            hXord = tmpX[:, hord]
        del(tmpX)
        return(hXord)

    def load_GTEx_eQTLs(self):
        self.qtlfile = self.qtldir + '/proc_qtls.cp.gz'
        try:
            # Try to load dataframe::
            print("[STATUS] Trying to load qtls from cp file")
            self.qtldf = AUX.load_pickle_gzip(self.qtlfile)
        except:
            print("[STATUS] Creating qtl dataframe")
            qtail = '.v7.signif_variant_gene_pairs.txt.gz'
            flist = glob.glob(self.qtldir + "/*" + qtail)
            self.qtldf = None
            for pairfile in flist:
                tissue = re.sub(qtail, "", re.sub(".*/", "", pairfile))
                print("- " + tissue)
                pdf = self.process_pairfile(pairfile)
                # Add tissue column:
                pdf['tissue'] = tissue
                if self.qtldf is None:
                    self.qtldf = pdf
                else:
                    self.qtldf = pd.concat([self.qtldf, pdf])
            print("[STATUS] Saving qtl dataframe")
            AUX.save_pickle_gzip(self.qtlfile, self.qtldf)
        self.assign_qtl_to_enh()

    def load_pchic(self):
        # Load PromoterCapture HiC links:
        self.lddir = self.datadir + '/linking_data/'
        self.pcfile = self.lddir + \
            'ActivePromoterEnhancerLinks_JavierrePChIC.tsv.gz'
        self.pcdf = pd.read_csv(self.pcfile, sep="\t")
        self.pcdf = self.pcdf.loc[:, ['baitID', 'oeID']].drop_duplicates()
        # Bait and overlap mappings:
        self.btdf = pd.read_csv(self.lddir + 'baitint.tsv', sep="\t")
        self.oedf = pd.read_csv(self.lddir + 'oeint.tsv', sep="\t")
        # Map these dfs to the clsids (bait is promoter):
        self.btdf = self.btdf.merge(self.tmdf)
        self.btdf = self.btdf.loc[:, ['baitID', 'cls']]
        self.btdf.columns = ['baitID', 'cls_promind']
        # Same for enhancers:
        self.oedf = self.oedf.merge(self.enhcoord)
        self.oedf = self.oedf.loc[:, ['oeID', 'cls']]
        self.oedf.columns = ['oeID', 'cls_enhind']
        # Merge back to matches:
        self.pcdf = self.pcdf.merge(self.btdf)
        self.pcdf = self.pcdf.merge(self.oedf)
        self.pcdf = self.pcdf.loc[:, ['cls_promind', 'cls_enhind']]
        self.pcdf = self.pcdf.drop_duplicates()
        print("[STATUS] Unique promoter - enhancer PCHiC interactions: " +
              fmtcm(len(self.pcdf)))
        # Merge with the test df:
        self.pcdf['ispc'] = 1
        self.testdf = self.testdf.merge(self.pcdf, how='left')
        self.testdf.ispc.fillna(0, inplace=True)
        print("[STATUS] Merge with test links: " +
              fmtcm(np.sum(self.testdf.ispc)) + " PCHiC links in " +
              fmtcm(len(self.testdf)) + " total possible links")

    def load_rmlinks(self):
        print("[STATUS] Loading Roadmap links")
        # Load Roadmap links:
        self.rdir = self.datadir + '/linking_data/RoadmapLinks/'
        self.rlfile = self.rdir + 'all_links_mapped.tsv.gz'
        self.rldf = pd.read_csv(self.rlfile, sep="\t")
        print("[STATUS] Start with: " + fmtcm(len(self.rldf)) + " links")
        # Same for enhancers:
        self.rldf = self.rldf.merge(self.enhcoord)
        self.rldf = self.rldf.loc[:, ['cls', 'GeneID']]
        self.rldf.columns = ['cls_enhind', 'short_gene']  # Diff gene version.
        self.tmdf['short_gene'] = [g.split('.')[0] for g in self.tmdf.GeneID]
        self.rldf = self.rldf.merge(self.tmdf)

        self.rldf = self.rldf.loc[:, ['cls_enhind', 'cls']]
        self.rldf.columns = ['cls_enhind', 'cls_promind']
        self.rldf.index = range(len(self.rldf))
        print("[STATUS] After merge, " + fmtcm(len(self.rldf)) + " links")
        # Merge with the test df:
        self.rldf['isrm'] = 1
        self.testdf = self.testdf.merge(self.rldf, how='left')
        self.testdf.isrm.fillna(0, inplace=True)
        print("[STATUS] Merge with test links: " +
              fmtcm(np.sum(self.testdf.isrm)) + " Roadmap links in " +
              fmtcm(len(self.testdf)) + " total possible links")

    # Process pairs of significant gene:
    def process_pairfile(self, pairfile):
        tmpdf = pd.read_csv(pairfile, sep="\t")
        # Get SNP location:
        fixcols = [n.split("_")[0:2] for n in tmpdf['variant_id'].tolist()]
        tmpdf['chr'] = [it[0] for it in fixcols]
        tmpdf['pos'] = [int(it[1]) for it in fixcols]
        # Return tissue, pairings
        keepcols = ['chr', 'pos', 'tss_distance', 'slope', 'gene_id']
        return(tmpdf[keepcols])

    def load_annotations(self):
        print("[STATUS] Loading annotations/metadata")
        # Load metadata - tissue to BSSID mapping
        self.anndf = pd.read_csv(
            self.anndir + '/updated_imputation_metadata.tsv', sep="\t")
        self.idmap = self.anndf[['BSSID', 'newgroup']]
        # TODO: Mapping of GTEx tissue to group tissue
        # Load all locations for DHS indexes:
        self.coordf = pd.read_csv(self.coordfile, sep="\t", header=None)
        self.coordf.columns = ['Chromosome', 'Start', 'End', 'name']
        # NOTE: coordf not in order of matrix - need coordnam:
        self.coordnam = pd.read_csv(self.coordnamfile, sep="\t")
        # Reorder coordf according to coordnam:
        print("[STATUS] Reordering coordinates list to match matrices")
        self.coordf = self.coordf.merge(self.coordnam)
        self.cord = np.argsort(self.coordf.cls, kind='mergesort')
        self.coordf = self.coordf.iloc[self.cord, :]
        self.coordf.cls = self.coordf.cls - 1  # Revert to 0-indexed
        # self.coordf.cls = range(len(self.coordf)) # TODO REMOVE
        self.coordf.index = range(len(self.coordf))
        self.coordf['mid'] = np.mean(self.coordf.loc[:, ['Start', 'End']],
                                     axis=1)
        self.midmap = self.coordf.loc[:, ['cls', 'mid']]
        # Add marginals:
        self.coordf['emarg'] = self.full_emarg
        self.coordf['pmarg'] = self.full_pmarg
        # Subset coordinates:
        self.enhcoord = self.coordf.iloc[self.enhind, :]
        self.promcoord = self.coordf.iloc[self.promind, :]
        self.enhcoord.index = range(len(self.enhcoord))
        self.promcoord.index = range(len(self.promcoord))
        # Turn into ranges object for intersections:
        self.enhgr = pr.PyRanges(self.enhcoord)
        self.promgr = pr.PyRanges(self.promcoord)
        # Load other annotations:
        self.load_gene_annotations()
        self.load_cluster_annotations()

    def load_gene_annotations(self, tol=0, load_gencode=False):
        print("[STATUS] Loading gene annotations")
        # GENEDF object from all GTEx tested genes:
        self.genedf = pd.read_csv(self.qtldir + "/all_egenes.txt", sep="\t")
        isrev = self.genedf.strand == '-'
        self.genedf['TSS'] = isrev * self.genedf.gene_end + \
            (-isrev) * self.genedf.gene_start
        # Make intersection object:
        self.genedf.gene_chr = ['chr' + n for n in self.genedf.gene_chr]
        self.tssdf = pd.DataFrame({'Chromosome': self.genedf.gene_chr,
                                   'Start': self.genedf.TSS - tol,
                                   'End': self.genedf.TSS + tol,
                                   'GeneID': self.genedf.gene_id,
                                   'Gene': self.genedf.gene_name})
        self.tssgr = pr.PyRanges(self.tssdf)
        if load_gencode:
            # GENCODEDF object from gencode annotation:
            print("[STATUS] Loading GENCODE annotations")
            self.gencodedf = pd.read_csv(
                self.anndir + '/gencode.gene.v30lift37.basic.annotation.gtf.gz',
                skiprows=5, header=None, sep="\t")
            self.gencodedf.columns = ['chr', 'group', 'type', 'start',
                                      'end', 'e1', 'strand', 'e2', 'info']
            isrev = self.gencodedf.strand == '-'
            self.gencodedf['TSS'] = isrev * self.gencodedf.end + \
                (-isrev) * self.gencodedf.start
            self.gcdf = pd.DataFrame({'Chromosome': self.gencodedf.chr,
                                      'Start': self.gencodedf.TSS - tol,
                                      'End': self.gencodedf.TSS + tol})
            self.gcgr = pr.PyRanges(self.gcdf)
        # TODO Re-organize dont put these functions here:
        # Assign promoters to genes (multiple prom per gene)
        self.assign_promoters_to_tss()
        # Create testdf:
        self.create_links_testset()

    def assign_promoters_to_tss(self, closest_cutoff=200):
        self.tmfile = self.qtldir + '/tss_map_closest.tsv'
        self.tmdf = pd.read_csv(self.tmfile, sep="\t", header=None)
        self.tmdf.columns = ['chr', 'tss', 'GeneID', 'name', 'dist_tm', 'tmid']
        # Remove all further than closest_cutoff:
        self.tmdf = self.tmdf.loc[self.tmdf.dist_tm <= closest_cutoff, :]
        # Remove all with no promoter marg:
        self.tmdf = self.tmdf.merge(self.coordf)
        self.tmdf = self.tmdf.loc[self.tmdf.pmarg > 0, :]
        self.genemap = self.tmdf.loc[:, ['GeneID', 'cls', 'tmid']]
        # Filtering only "promoter" identified:
        self.ptmdf = self.genemap.merge(self.promcoord)
        self.genemap.columns = ['GeneID', 'cls_promind', 'tmid']
        # NOTE: 19k vs 15k promoters

    # Create test links:
    def create_links_testset(self):
        print("[STATUS] Making dataframe for testing E-P pairs")
        # Promoter regions +/- 500kb
        self.pfile = self.qtldir + '/all_1MB_intersections.tsv.gz'
        self.pairsdf = pd.read_csv(self.pfile, sep="\t", header=None)
        self.pairsdf.columns = ['name', 'tmid']
        # GeneID to promoter index:
        self.pairsdf = self.pairsdf.merge(self.genemap)
        # NOTE: Keep only enhancer chunks:
        self.pairsdf = self.pairsdf.merge(self.enhcoord)
        self.pairsdf = self.pairsdf.loc[:, ['cls_promind', 'cls']
                                        ].drop_duplicates()
        self.pairsdf.columns = ['cls_promind', 'cls_enhind']
        print("[STATUS] Created links testset with " +
              fmtcm(len(self.pairsdf)) + " possible links in +/- 1MB")

    def assign_qtl_to_enh(self):
        print("[STATUS] Assigning qtls to enhancers")
        self.qtldf.index = range(len(self.qtldf))
        self.qtldf.chr = ['chr' + n for n in self.qtldf.chr]
        # Load the intersection assignments:
        self.qifile = self.qtldir + '/int_snplist_ann.tsv'
        self.qidf = pd.read_csv(self.qifile, sep="\t", header=None)
        self.qidf.columns = ['chr', 'pos', 'dist', 'cls', 'name']
        # Filter by some cutoff distance:
        cutoff_dist = 1000
        self.qidf = self.qidf.loc[self.qidf.dist <= cutoff_dist, :]
        self.qidf.index = range(len(self.qidf))
        self.qtldf = self.qtldf.merge(self.qidf)
        # From chunk name and snp position get QTL links:
        self.closeqtl = self.qtldf.loc[:, ['name', 'gene_id', 'cls', 'tissue']]
        self.closeqtl = self.closeqtl.merge(self.enhcoord)
        self.closeqtl = self.closeqtl.loc[:, ['name', 'gene_id',
                                              'cls', 'tissue']]
        self.closeqtl.columns = ['name', 'GeneID', 'cls_enhind', 'tissue']
        self.closeqtl = self.closeqtl.merge(self.genemap)
        self.qtltisdf = self.closeqtl.loc[:, ['cls_promind',
                                              'cls_enhind', 'tissue']].copy()
        self.closeqtl = self.closeqtl.loc[:, ['cls_promind', 'cls_enhind']]
        self.closeqtl = self.closeqtl.drop_duplicates()
        print("[STATUS] Unique promoter - enhancer reduced QTLs to " +
              fmtcm(len(self.closeqtl)) + " QTL links")
        # Merge with the test df:
        self.closeqtl['isqtl'] = 1
        self.testdf = self.pairsdf.merge(self.closeqtl, how='left')
        self.testdf.isqtl.fillna(0, inplace=True)
        print("[STATUS] Merge with test links creates final dataset: " +
              fmtcm(np.sum(self.testdf.isqtl)) + " QTL links in " +
              fmtcm(len(self.testdf)) + " total possible links")

    def load_cluster_annotations(self):
        # Load cluster assignments:
        ctag = 'cls_merge2_wH3K27ac100_' + str(self.nclust) + '_assignments'
        # NOTE: id column is 1-indexed in the bedfile:
        self.eclsdf = pd.read_csv(
            self.datadir + '/' + ctag + '_ENH_052819.bed', sep="\t")
        self.pclsdf = pd.read_csv(
            self.datadir + '/' + ctag + '_PROM_052819.bed', sep="\t")
        self.eclsdf.columns = ['id_enh', 'cls_enhind']
        self.pclsdf.columns = ['id_prom', 'cls_promind']
        # Revert to 0-indexed
        self.eclsdf.cls_enhind = self.eclsdf.cls_enhind - 1
        self.pclsdf.cls_promind = self.pclsdf.cls_promind - 1

    def match_clusters(self, cutoff=0.25):
        # Get cluster centers for enhancers and promoters:
        emat = self.eX.multiply(1 * (self.h1X > 0))[self.enhind, :]
        pmat = self.pX.multiply(1 * (self.h1X > 0))[self.promind, :]
        self.ecent = self.calculate_centers(emat, self.eclsdf.id_enh)
        self.pcent = self.calculate_centers(pmat, self.pclsdf.id_prom)
        # Calculate distance matrix between cluster centers:
        Zec = 1 * (self.ecent > cutoff)
        Zpc = 1 * (self.pcent > cutoff)
        self.simep, self.intcep = matrix_jaccard(Zec, Zpc)
        mat = coo_matrix(self.simep)
        self.moddf = pd.DataFrame({'id_enh': mat.row,
                                   'id_prom': mat.col,
                                   'mod_sim': 1 - mat.data})
        # self.reord = AUX.reorder_by_linkage(self.simep)

    def calculate_centers(self, mat, assign):
        ua, c = np.unique(assign, return_counts=True)
        # print(ua)
        centers = np.zeros((mat.shape[1], len(ua)))
        for i in ua:
            slc = np.where(assign == i)[0]
            centers[:, i] = np.array(np.mean(mat[slc, :], 0))[0]
        return(centers)


# Tools for overlap:
def overlap_bothix(scdf, ocdf):
    it = ncls.NCLS(
        ocdf.Start.values, ocdf.End.values, ocdf.index.values)
    ind1, ind2 = it.all_overlaps_both(
        scdf.Start.values, scdf.End.values, scdf.index.values)
    return(ind1, ind2)


def reindex(gr):
    df = gr.df
    df.index = range(len(df))
    return(pr.PyRanges(df))


def collate_idlists(idlists, col1, col2):
    testdf = None
    for chrom in idlists.keys():
        tdf = pd.DataFrame({col1: idlists[chrom][0],
                            col2: idlists[chrom][1]})
        if testdf is None:
            testdf = tdf
        else:
            testdf = testdf.append(tdf)
    testdf.index = range(len(testdf))
    return(testdf)


# Given X and Y that have the same rows (epigenomes):
def matrix_jaccard(X, Y):
    intrsct = X.T.dot(Y * 1.0)
    csX = np.squeeze(np.array(np.sum(X, axis=0)))
    csY = np.squeeze(np.array(np.sum(Y, axis=0)))
    unions = csX[:, None] + csY - intrsct
    dist = 1.0 - intrsct / unions
    return(dist, intrsct)


def fmtcm(x):
    return("{:,}".format(int(x)))


def get_data(pref):
    with gzip.open(pref + "_csr.cp.gz", 'rb') as f:
        X = pickle.load(f, encoding='latin1')
    return(X)


def get_names(pref):
    with gzip.open(pref + "_attr.cp.gz", 'rb') as f:
        attr = pickle.load(f, encoding='latin1')
    return(attr['names'])


if __name__ == "__main__":
    fire.Fire(load_linking_data)
