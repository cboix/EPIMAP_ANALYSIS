#!/bin/bash
# ---------------------------------------
# Reduce the motif matches to BKG + DHSs:
# ---------------------------------------
start=`date +%s`
hostname -f

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

export NSHUFFLES=100
export K=8

if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [OPTIONS] 
    -d     Dataset to use (defaults to collated JASPAR/HOCOMOCO/Jolma2013 motifs)
    -k     KMER length (for p-value cutoff) [default: $K]
    -s     Number of shuffles to calculate [default: $NSHUFFLES]"
    exit 1
fi

while getopts d:k:s: o
do      case "$o" in
    d)		export ALL_PFM="$OPTARG";;
    k)		export K="$OPTARG";;
    s)      export NSHUFFLES="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done


# Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
MOTIFDIR=$DBDIR/motifs
MATCHDIR=$MOTIFDIR/matches
mkdir -p $MOTIFDIR $MATCHDIR

cd $MOTIFDIR

export TASK=${SGE_TASK_ID}
echo "[STATUS] TASK NUMBER $TASK"
TMP_DIR=${TMP}/motif_scans_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR}

# Genome to use:
GENOME=$MAINDIR/genomes/hg19/hg19.fa

# Motifs dataset:
if [[ -z "$ALL_PFM" ]]; then
    ALL_PFM=$MOTIFDIR/collated_pfm.txt
    ALL_MOT=$MOTIFDIR/collated_motifs.txt
    echo "[STATUS] Defaulting to dataset: ${ALL_PFM}"
fi

# --------------
# Extract motif:
# --------------
awk -v num=$TASK '$1 ~ />/{a=a+1}{if(a == num){print $0}}' $ALL_PFM > ${TMP_DIR}/pfm.txt
MOTNAM=$( awk '$1 ~ />/{sub(">","",$1); print $1}' ${TMP_DIR}/pfm.txt )
echo "Motif is: $MOTNAM"

# Get motif:
awk -v mot=$MOTNAM '$1 ~ />/{a = 0}$1 == ">"mot{a=1}{if(a == 1){print $0}}' $ALL_MOT > ${TMP_DIR}/mot.txt


KMERS=8
KEEPWITHMAX=0
PATH="$PATH:$MAINDIR/MOTIF_VALIDATION/bin"
SCORE=$( MotifAddPValCutoff.sh -k "$KMERS" -s "_%dmer" -m $KEEPWITHMAX $TMP_DIR/mot.txt | awk 'NR==1{print $2}' )
echo "TFM-pvalue score is $SCORE"

# Reduce the matches to BKG and DHS regions (CORECOORD):
BACKREG=${MOTIFDIR}/back-regions.txt.gz
MOTPREF=${MATCHDIR}/matches_wshuffles_n${NSHUFFLES}_${MOTNAM}
FINALFILE=${MOTPREF}.tsv.gz

# Merged regions:
BEREG=${MOTIFDIR}/back-enh-regions.txt

# ----------------------------
# Notes on background regions:
# ----------------------------
# Length of core coords:
# awk -vOFS="\t" '{a=$3 - $2 + a}END{print a}' $CORECOORD
# 728977213
# Check how much of corecoords is contained in back regions:
# bedtools intersect -a ${BACKREG} -b ${CORECOORD} | awk -vOFS="\t" '{a=$3 - $2 + a}END{print a}'
# 440450075
# bedtools intersect -a ${CORECOORD} -b ${BACKREG} | awk -vOFS="\t" '{a=$3 - $2 + a}END{print a}'
# 440479311
# Intergenic regions aren't the best control, because 2.4M elements do match but 838k dont at all (and 250k+ match 2+ times) - genic DHS/enhancers.
# bedtools intersect -a ${CORECOORD} -b ${BACKREG} -c | cut -f5 | sort | uniq -c

# Background comparisons to evaluate:
# 1. Intergenic + DHS sites (feed in bkg counts to enricher)
# 2. Full genome (feed in counts directly to enricher)
# 3. (TODO?) Crop regions to Intergenic DHS only.... --> Leads to issues down the road.

# -------------------------------------------
# Create the intergenic + DHS sites bed file:
# -------------------------------------------
if [[ ! -s $BEREG.gz ]]; then
    zcat $BACKREG > ${BEREG}
    awk -vOFS="\t" '{print $1,$2,$3}' $CORECOORD >> ${BEREG}
    sort -k1,1V -k2,2n ${BEREG} | bedtools merge -i - | gzip -c > ${BEREG}.gz
    # BKGLEN=$( zcat ${BEREG}.gz | awk -vOFS="\t" '{a=$3 - $2 + a}END{print a}' )
    # 1773347194
fi
BKGLEN="1773347194"

# Outputs:
# REDFILE=${MOTPREF}.reduced.tsv.gz
DHSFILE=${MOTPREF}.onlyDHS.tsv.gz
BKGFILE=${MOTPREF}.bkg.tsv.gz
BEFILE=${MOTPREF}.bkgdhs.tsv.gz # Joint background
ERESFILE=${MOTPREF}.bkgdhs.epi.enrich.tsv.gz
CRESFILE=${MOTPREF}.bkgdhs.cls.enrich.tsv.gz
# Subsetted by epi/screen/both
CDRESFILE=${MOTPREF}.bkgdhs.screen_dELS.cls.enrich.tsv.gz
CNRESFILE=${MOTPREF}.bkgdhs.screen_non_dELS.cls.enrich.tsv.gz
CERESFILE=${MOTPREF}.bkgdhs.epimap_non_screen.cls.enrich.tsv.gz
# In epigenomes:
EDRESFILE=${MOTPREF}.bkgdhs.screen_dELS.epi.enrich.tsv.gz
ENRESFILE=${MOTPREF}.bkgdhs.screen_non_dELS.epi.enrich.tsv.gz
EERESFILE=${MOTPREF}.bkgdhs.epimap_non_screen.epi.enrich.tsv.gz

if [[ ! -s $ERESFILE ]] || [[ ! -s $CRESFILE ]] || [[ ! -s $EDRESFILE ]]; then
    # Turn matches into bedfile for bedtools intersect
    # zcat $FINALFILE | awk -vOFS="\t" -v score=$SCORE '$5 > score{print $1,$2,$2,$4,NR}' > ${TMP_DIR}/matches.bed
    zcat $FINALFILE | awk -vOFS="\t" '{print $1,$2,$2,$4,NR}' > ${TMP_DIR}/matches.bed
    # Get counts overall
    cut -f4 ${TMP_DIR}/matches.bed | sort | uniq -c > ${TMP_DIR}/pfms_full_counts.txt

    # Add the number of bkg regions too:
    # TODO: Count bkg LENGTH, not regions (div by avg. DHS length)?

    # Get matches in the DHSs only and in background region sets overall:
    if [[ ! -s $DHSFILE ]] || [[ $FINALFILE -nt $DHSFILE ]]; then 
        bedtools intersect -a ${TMP_DIR}/matches.bed -b ${CORECOORD} -wa -wb | awk -vOFS="\t" '{print $4,$5,$9}' | gzip -c > ${DHSFILE}
        zcat $DHSFILE | wc
    fi

    if [[ ! -s $BKGFILE ]] || [[ $FINALFILE -nt $BKGFILE ]]; then
        bedtools intersect -a ${TMP_DIR}/matches.bed -b ${BACKREG} -wa -wb | awk -vOFS="\t" '{print $4,$5}' | gzip -c > ${BKGFILE}
        zcat $BKGFILE | wc
    fi

    if [[ ! -s $BEFILE ]] || [[ $FINALFILE -nt $BEFILE ]]; then
        bedtools intersect -a ${TMP_DIR}/matches.bed -b ${BEREG} -wa -wb | awk -vOFS="\t" '{print $4,$5}' | gzip -c > ${BEFILE}
        zcat $BEFILE | wc
    fi

    # Count number of matches in DHSs and background regions.
    zcat ${DHSFILE} | cut -f1 | sort | uniq -c > ${TMP_DIR}/pfms_dhs_counts.txt
    zcat ${BKGFILE} | cut -f1 | sort | uniq -c > ${TMP_DIR}/pfms_bkg_counts.txt
    zcat ${BEFILE} | cut -f1 | sort | uniq -c > ${TMP_DIR}/pfms_bkgdhs_counts.txt

    # Make table comparing fractional enrichments:
    join -1 2 -2 2 ${TMP_DIR}/pfms_bkgdhs_counts.txt ${TMP_DIR}/pfms_full_counts.txt | join -1 2 -2 1 ${TMP_DIR}/pfms_bkg_counts.txt - | join -1 2 -2 1 ${TMP_DIR}/pfms_dhs_counts.txt - | awk -vOFS="\t" 'BEGIN{print "motif","ndhs","nbkg","nbkgdhs","nfull", "dhsvbkg","dhsvbkgdhs","dhsvfull"}{print $1,$2,$3,$4,$5,$2/$3,$2/$4,$2/$5}' | sort -n -k7 > ${MOTPREF}.DHSfracs.txt
    cat ${MOTPREF}.DHSfracs.txt

    # ----------------------------------------------------------
    # Quick enrichments against the epigenomes and the clusters:
    # Using all DHS sites as a background
    # ----------------------------------------------------------
    # Mapping files for epigenomes and clusters:
    MODELDIR=$DBDIR/ChromHMM/calls/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust
    EPIFILE=${MODELDIR}/cls_merge2_wH3K27ac100_raw/motif_enrichment/observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_full.bed.gz
    CLSFILE=${MODELDIR}/cls_merge2_wH3K27ac100_300/motif_enrichment/observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_full.bed.gz

    # Calculate the fraction motif enrichment in/out, given DHS file, matches files, and mapping file. Uses the all DHS hits as background
    if [[ ! -s $CRESFILE ]]; then
        python $BINDIR/motif_enrichment/compute_motif_enrichment.py main --regfile ${DHSFILE} --mapfile ${CLSFILE} --outprefix ${MOTPREF}.bkgdhs.cls --motif $MOTNAM --bkgcountfile ${MOTPREF}.DHSfracs.txt --bkglen ${BKGLEN}
        # python $BINDIR/motif_enrichment/compute_motif_enrichment.py main --regfile ${DHSFILE} --mapfile ${CLSFILE} --outprefix ${MOTPREF}.cls --motif $MOTNAM 
    fi

    # --------------------------------------------------
    # Add specific intersections by diff. enhancer sets:
    # --------------------------------------------------
    # TODO: If time, need to do intersections for "screen only" (very small set)
    # CDRESFILE=${MOTPREF}.bkgdhs.screen_dELS.cls.enrich.tsv.gz
    if [[ ! -s $CDRESFILE ]]; then
        echo "[STATUS] Running against intersection of dELS + Epimap enhancers"
        python $BINDIR/motif_enrichment/compute_motif_enrichment.py main --regfile ${DHSFILE} --mapfile ${CLSFILE} --outprefix ${MOTPREF}.bkgdhs.screen_dELS.cls --motif $MOTNAM --bkgcountfile ${MOTPREF}.DHSfracs.txt --bkglen ${BKGLEN} --elemlist $DBDIR/epimap_screen_dELS.txt
    fi

    # CNRESFILE=${MOTPREF}.bkgdhs.screen_non_dELS.cls.enrich.tsv.gz
    # if [[ ! -s $CNRESFILE ]]; then
    #     echo "[STATUS] Running against intersection of NON-dELS + Epimap enhancers"
    #     python $BINDIR/motif_enrichment/compute_motif_enrichment.py main --regfile ${DHSFILE} --mapfile ${CLSFILE} --outprefix ${MOTPREF}.bkgdhs.screen_non_dELS.cls --motif $MOTNAM --bkgcountfile ${MOTPREF}.DHSfracs.txt --bkglen ${BKGLEN} --elemlist $DBDIR/epimap_screen_non_dELS.txt
    # fi

#     # CERESFILE=${MOTPREF}.bkgdhs.epimap_non_screen.cls.enrich.tsv.gz
#     if [[ ! -s $CERESFILE ]]; then
#         echo "[STATUS] Running against unique Epimap enhancers (non-SCREEN)"
#         python $BINDIR/motif_enrichment/compute_motif_enrichment.py main --regfile ${DHSFILE} --mapfile ${CLSFILE} --outprefix ${MOTPREF}.bkgdhs.epimap_non_screen.cls --motif $MOTNAM --bkgcountfile ${MOTPREF}.DHSfracs.txt --bkglen ${BKGLEN} --elemlist $DBDIR/epimap_non_screen_elements.txt
#     fi

    # ---------------------------
    # Epigenomes + intersections:
    # ---------------------------
    if [[ ! -s $ERESFILE ]]; then
        python $BINDIR/motif_enrichment/compute_motif_enrichment.py main --regfile ${DHSFILE} --mapfile ${EPIFILE} --outprefix ${MOTPREF}.bkgdhs.epi --motif $MOTNAM --bkgcountfile ${MOTPREF}.DHSfracs.txt --bkglen ${BKGLEN} --bkgcountfile ${MOTPREF}.DHSfracs.txt --bkglen ${BKGLEN}
    fi

    if [[ ! -s $EDRESFILE ]]; then
        echo "[STATUS] Running against intersection of dELS + Epimap enhancers"
        python $BINDIR/motif_enrichment/compute_motif_enrichment.py main --regfile ${DHSFILE} --mapfile ${EPIFILE} --outprefix ${MOTPREF}.bkgdhs.screen_dELS.epi --motif $MOTNAM --bkgcountfile ${MOTPREF}.DHSfracs.txt --bkglen ${BKGLEN} --elemlist $DBDIR/epimap_screen_dELS.txt
    fi

    # if [[ ! -s $ENRESFILE ]]; then
    #     echo "[STATUS] Running against intersection of NON-dELS + Epimap enhancers"
    #     python $BINDIR/motif_enrichment/compute_motif_enrichment.py main --regfile ${DHSFILE} --mapfile ${EPIFILE} --outprefix ${MOTPREF}.bkgdhs.screen_non_dELS.epi --motif $MOTNAM --bkgcountfile ${MOTPREF}.DHSfracs.txt --bkglen ${BKGLEN} --elemlist $DBDIR/epimap_screen_non_dELS.txt
    # fi

    # if [[ ! -s $EERESFILE ]]; then
    #     echo "[STATUS] Running against unique Epimap enhancers (non-SCREEN)"
    #     python $BINDIR/motif_enrichment/compute_motif_enrichment.py main --regfile ${DHSFILE} --mapfile ${EPIFILE} --outprefix ${MOTPREF}.bkgdhs.epimap_non_screen.epi --motif $MOTNAM --bkgcountfile ${MOTPREF}.DHSfracs.txt --bkglen ${BKGLEN} --elemlist $DBDIR/epimap_non_screen_elements.txt
    # fi


else
    echo "[STATUS] File already exists at: $ERESFILE and $CRESFILE"
    ls -sh $ERESFILE $CRESFILE
fi


rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished motif matches ($TASK) sucessfully in $runtime seconds."
