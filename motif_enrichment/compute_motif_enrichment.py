#!/usr/bin/python
# -------------------------------------------------
# Compute motif enrichment.
# From region overlaps + mapping to sets of regions
# -------------------------------------------------
import os
import re
import fire
import time
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

rng = np.random


class compute_motif_enrichment(object):
    def __init__(self, regfile, mapfile, outprefix, elemlist=None,
                 bkgcountfile=None, bkglen=3.2e9, bkgstr='nbkgdhs',
                 mergedups=False, confidz=1.5, motif=None):
        # Arguments:
        self.regfile = regfile
        self.mapfile = mapfile
        self.outprefix = outprefix
        self.outfile = self.outprefix + ".enrich.tsv.gz"
        self.elemlist = elemlist
        # Background parameters, case where bkgcounts provided:
        self.bkgcountfile = bkgcountfile
        self.bkgstr = bkgstr
        self.bkglen = bkglen
        # Other parameters:
        self.mergedups = mergedups # Merge same hit to two regions or not
        self.confidz = confidz
        self.motif = motif

    def main(self):
        self.read_data()
        if self.bkgcountfile is not None:
            self.read_bkg()
        self.calculate_duplicate_correction()
        self.merge_with_mapping()
        self.count_hits_motif()
        self.calculate_enrichments()
        self.write_enrichments()

    def read_data(self):
        t1 = time.time()
        # Motif matches:
        print("[STATUS] Reading matches to regions file")
        self.rgdf = pd.read_csv(self.regfile, skiprows=0, header=None,
                                sep="\t", names=['motif','matchnum','chunk'])
        print(self.rgdf.shape)
        if self.elemlist is not None:
            print("[STATUS] Reading element subset list")
            self.elems = pd.read_csv(self.elemlist,
                                     skiprows=0, header=None).to_numpy().T[0]
            print(len(self.elems), "elements read")
            print("[STATUS] Reducing to subset of elements")
            self.rgdf = self.rgdf.loc[self.rgdf['chunk'].isin(self.elems)]
            print(self.rgdf.shape)
        # Mapping file:
        print("[STATUS] Reading mapping file")
        self.mapdf = pd.read_csv(self.mapfile, skiprows=0, header=None,
                                 sep="\t", names=['chunk','set'])
        print(self.mapdf.shape)
        if self.elemlist is not None:
            print("[STATUS] Reducing to subset of elements")
            self.mapdf = self.mapdf.loc[self.mapdf['chunk'].isin(self.elems)]
            print(self.mapdf.shape)
        print("[STATUS] Finished reading inputs in " +
              str(round(time.time() - t1)) + 's')
        # Count number of foreground, background regions for mapfile:
        print("[STATUS] Calculating the raw region counts")
        self.mapdf['regtot'] = 1
        self.mapmarg = self.mapdf.groupby(
            'set', as_index=False)['regtot'].agg('sum')
        # Approximate number of DHS-equivalent regions
        if self.bkgcountfile is None:
            self.maptotal = np.sum(self.mapmarg['regtot'])
            # NOTE: Overlapping counting if epigenomes, shouldn't do this.
        else:
            self.maptotal = np.floor(self.bkglen / 200)
        print("[STATUS] Raw number of background regions:", self.maptotal)

    def read_bkg(self):
        print("[STATUS] Reading background counts file")
        self.bgdf = pd.read_csv(self.bkgcountfile, skiprows=0, sep="\t")
        print(self.bgdf.head())
        self.bgdf = self.bgdf[['motif', self.bkgstr]]
        self.bgdf.columns = ['motif', 'total']
        self.bgdf['is_motif'] = (self.bgdf.motif == 'pfm')
        # Create totals from these background counts:
        self.totdf = self.bgdf.groupby(
            ['is_motif'], as_index=False)['total'].agg('sum')

    def calculate_duplicate_correction(self):
        if self.mergedups:
            # Get the margin for the rgdf data and merge back to rgdf:
            print("[STATUS] Calculating the duplicate corrections")
            self.rgdf['pcount'] = 1
            self.rgmarg = self.rgdf.groupby(
                'matchnum', as_index=False)['pcount'].agg('sum')
            self.rgmarg.columns = ['matchnum','total']
            self.rgmarg['weight'] = 1 / self.rgmarg['total']
            self.rgdf = self.rgdf.merge(self.rgmarg)
            print(self.rgdf.shape)
        else:
            self.rgdf['weight'] = 1

    def merge_with_mapping(self):
        print("[STATUS] Merging hits with region mapping")
        self.rgdf_short = self.rgdf[['motif','chunk','weight']]
        self.rgdf_short = self.rgdf_short.merge(self.mapdf)

    def count_hits_motif(self):
        # Aggregate number of hits per motif set:
        print("[STATUS] Counting hits by motif:")
        self.rgdf_short['is_motif'] = (self.rgdf_short.motif == 'pfm')
        self.ctdf = self.rgdf_short.groupby(
            ['is_motif','set'], as_index=False)['weight'].agg('sum')
        # Overall counts:
        if self.bkgcountfile is None:
            print("[STATUS] Using internal background counts")
            self.totdf = self.ctdf.groupby(
                ['is_motif'], as_index=False)['weight'].agg('sum')
            self.totdf.columns = ['is_motif','total']
        else:
            print("[STATUS] Using provided background counts")
        self.ctreddf = self.ctdf.merge(self.totdf)
        self.pfmdf = self.ctreddf[self.ctreddf.is_motif == True]
        self.shfdf = self.ctreddf[self.ctreddf.is_motif == False]
        # Print outputs:
        self.pfmdf.columns = ['pfmcol', 'set','pfg','pbg']
        self.shfdf.columns = ['shfcol', 'set','cfg','cbg']
        self.pfmdf = self.pfmdf[['set','pfg','pbg']].merge(self.shfdf[['set','cfg','cbg']], how='outer')
        self.pfmdf.fillna(value={'pfg': 0, 'pbg': self.pfmdf.pbg[0]})
        # Round down all foreground; up all background:
        self.pfmdf.pfg = self.pfmdf['pfg'].apply(np.floor)
        self.pfmdf.cfg = self.pfmdf['cfg'].apply(np.floor)
        self.pfmdf.pbg = self.pfmdf['pbg'].apply(np.ceil)
        self.pfmdf.cbg = self.pfmdf['cbg'].apply(np.ceil)
        self.pfmdf = self.pfmdf.merge(self.mapmarg)
        print(self.pfmdf.shape)

    def calculate_enrichments(self):
        # M pop, n tot, N draws, x hits
        self.pfmdf['pval_raw'] = -np.log10(hypergeom.sf(
            # Orig, fixed:
            # self.pfmdf.pfg - 1, self.maptotal, self.pfmdf.cfg, self.pfmdf.regtot))
            self.pfmdf.pfg - 1, self.maptotal, self.pfmdf.pbg, self.pfmdf.regtot))
        self.pfmdf['pval_ctrl'] = -np.log10(hypergeom.sf(
            self.pfmdf.pfg - 1, self.pfmdf.cbg, self.pfmdf.cfg, self.pfmdf.pbg))
        # TODO must be flipped
        self.pfrac = self.calc_binomial_CI(
            self.pfmdf.pfg, self.pfmdf.pbg, self.confidz, max=False)
        self.cfrac = self.calc_binomial_CI(
            self.pfmdf.cfg, self.pfmdf.cbg, self.confidz, max=True)
        self.pfrac = self.pfrac.fillna(0)
        self.cfrac = self.cfrac.fillna(0)
        self.pfrac[self.pfrac < 0] = 0
        self.cfrac[self.cfrac < 0] = 0
        self.rctrl = np.log2(self.pfrac / self.cfrac)
        self.rctrl[self.pfrac == 0] = 0
        # TODO: What if highly depleted: should switch min and max for CI.
        self.pfmdf['pfrac'] = self.pfrac
        self.pfmdf['cfrac'] = self.cfrac
        self.pfmdf['ratio_ctrl'] = self.rctrl
        print('Max enrichment (log2)', round(np.max(self.pfmdf['ratio_ctrl']),3))
        print('Max p-value raw (-log10)', round(np.max(self.pfmdf['pval_raw']),3))
        print('Max p-value ctrl (-log10)', round(np.max(self.pfmdf['pval_ctrl']),3))
        if self.motif is not None:
            self.pfmdf['motif'] = self.motif

    def write_enrichments(self):
        self.pfmdf.to_csv(self.outfile, sep="\t", index=False, float_format="%0.4e")

    # Wilson_score_interval from:
    # https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    def calc_binomial_CI(self, ns, n, z, max=True):
        nf = n - ns
        zsq = z**2
        p = (ns + zsq/2) / (n + zsq)
        pm = (z / (n + zsq)) * np.sqrt((ns * nf) / n + zsq / 4)
        return(p + (2 * max - 1) * pm)


if __name__ == "__main__":
    fire.Fire(compute_motif_enrichment)
