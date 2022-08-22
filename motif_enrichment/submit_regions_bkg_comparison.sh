#!/bin/bash
# ---------------------------------
# Submit motif enrichment pipeline:
# With example for scATAC here
# ---------------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

# Variables and files:
MOTIFDIR=$DBDIR/motifs
ALL_PFM=$MOTIFDIR/collated_pfm.txt
MNUM=$( grep ">" $ALL_PFM | wc -l )

conda activate base

# scATAC example run:
clist=( Ast Excitatory Inhibitory Microglia OPC Oligo Per_Endo )

ATDIR=${MAINDIR}/DEVTRAJ/db/AD_scATAC
DDIR=${ATDIR}/DiffPeaks100520

for cell in ${clist[@]}; do 
    echo "$cell"
    adfile=${DDIR}/${cell}.binomial.ADbyBraak.Diff.sort.txt
    ctfile=${DDIR}/CelltypeSpecif.peak.${cell}.fixed.txt
    outfile=${DDIR}/AD.vs.CT.counts.${cell}.tsv
    # ovlfile=${DDIR}/AllCTpeaks.txt # Alternative to Index DHSs
    ovlfile=${CORECOORD}
    out2file=${DDIR}/CT.vs.DHS.counts.${cell}.tsv
    wc -l $adfile
    wc -l $ctfile

    # Difftl enrichments against cell-type specific:
    qsub -cwd -l h_vmem=20G -l h_rt=8:00:00 -N motif_counts_difftl_${cell} -j y -b y -V -r y -o $DBDIR/out/ "$BINDIR/motif_enrichment/run_regions_bkg_comparison.sh -i ${adfile} -o ${outfile} -b ${ctfile}"

    # Cell-type specific enrichments against DHSs
    qsub -cwd -l h_vmem=20G -l h_rt=8:00:00 -N motif_counts_ovl_${cell} -j y -b y -V -r y -o $DBDIR/out/ "$BINDIR/motif_enrichment/run_regions_bkg_comparison.sh -i ${ctfile} -o ${out2file} -b ${ovlfile}"
done

