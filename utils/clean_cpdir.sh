#!/bin/bash
# For cleaning out the preprocessed and aggregated matrices
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES

CDIR=$CALLDIR/observed_aux_18_on_mixed_impobs_QCUT/

rm ${TMP}/dir.tsv
echo $CHMM_FMTDIR/DNase-seq/cpfiles >> $TMP/dir.tsv
echo $CHMM_FMTDIR/H3K27ac/cpfiles >> $TMP/dir.tsv
echo $CHMM_FMTDIR/H3K4me1/cpfiles >> $TMP/dir.tsv
echo $CHMM_FMTDIR/H3K4me2/cpfiles >> $TMP/dir.tsv
echo $CHMM_FMTDIR/H3K4me3/cpfiles >> $TMP/dir.tsv
echo $CHMM_FMTDIR/H3K9ac/cpfiles >> $TMP/dir.tsv
echo $CDIR/STATEBYLINE/cpfiles >> $TMP/dir.tsv

# TODO: CHECK THE FINAL OUTPUTS

while read cpdir; do 
    cd $cpdir
    echo $cpdir
    while read chrom size; do 
        echo $chrom
        rm BSS*_${chrom}_*nonovl*csr.cp.gz
        rm BSS*_${chrom}_*nonovl*attr.cp.gz
        rm BSS*_${chrom}_*csr.cp.gz
        rm BSS*_${chrom}_*attr.cp.gz
    done < $CHROMSIZES
done < $TMP/dir.tsv

# Remove the raw signal cp files:
cd $IMPUTED_DIR/cpfiles;
while read chrom size; do 
    echo $chrom
    while read mark; do 
        echo $mark
        # FOR IMPDIR
        # rm ${chrom}*_BSS*_${mark}*nonovl*csr.cp.gz
        # rm ${chrom}*_BSS*_${mark}*nonovl*attr.cp.gz
        rm ${chrom}*_BSS*_${mark}*csr.cp.gz
        rm ${chrom}*_BSS*_${mark}*attr.cp.gz
    done < $MARKS_LIST
done < $CHROMSIZES

cd $CONVERTED_DIR/cpfiles;
while read chrom size; do 
    echo $chrom
    while read mark; do 
        echo $mark
        # FOR IMPDIR
        # rm ${chrom}*_${mark}_BSS*nonovl*csr.cp.gz
        # rm ${chrom}*_${mark}_BSS*nonovl*attr.cp.gz
        rm ${chrom}*_${mark}_BSS*csr.cp.gz
        rm ${chrom}*_${mark}_BSS*attr.cp.gz
    done < $MARKS_LIST
done < $CHROMSIZES

# Clean out compiled tables
while read mark; do 
    cpdir=$CONVERTED_DIR/$mark
    if [[ -d $cpdir ]]; then
        echo $mark
        cd $cpdir
        rm *.cp.gz
        rm *.npz
    fi
    cpdir=$IMPUTED_DIR/$mark
    if [[ -d $cpdir ]]; then
        echo $mark
        cd $cpdir
        rm *.cp.gz
        rm *.npz
    fi
    cpdir=$CHMM_FMTDIR/$mark
    if [[ -d $cpdir ]]; then
        echo $mark
        cd $cpdir
        rm *.cp.gz
        rm *.npz
    fi
done < $MARKS_LIST

cd $CDIR/STATEBYLINE/ENH
rm *.cp.gz
rm *.npz

cd $CDIR/STATEBYLINE/PROM
rm *.cp.gz
rm *.npz

