#!/bin/bash
# -----------------------------------------------------
# 1. Extract mark (H3K27ac) around peaks (of DNase-seq)
#    for each epigenome
# 2. Plot each epigenome
# 3. Collect aggregate statistics
# 4. Plot aggregate statistics against 
#    stats from masterlist calls
# -----------------------------------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

TMP_DIR=${TMP}/extract_bw_peaks_general_${RANDOM}
mkdir -p ${TMP_DIR}

# NOTE: Using the mixobs table to pick pairs for each sample:
# TODO: Allow both, but give diff. resultdir?
MAININFO=$MIXOBS_TAB
PKMARK="DNase-seq"
BWMARK="H3K27ac"
EXTEND=250
MDPT=6

RESULTDIR=$PKIMP_DIR/${BWMARK}_inpeaksof_${PKMARK}/
mkdir -p ${RESULTDIR}

# -----------------
# 0. Make infofile:
# -----------------
INFOFILE=${TMP_DIR}/infofile.tsv
grep "${PKMARK}" ${MAININFO} | awk -vOFS="\t" '{print $1,$3"_peaks.summit.narrowPeak.gz"}' | sort +0 -1 > ${TMP_DIR}/peaks.tsv
grep "${BWMARK}" ${MAININFO} | awk -vOFS="\t" '{print $1,$3}' | sort +0 -1  > ${TMP_DIR}/bwpref.tsv
join ${TMP_DIR}/peaks.tsv ${TMP_DIR}/bwpref.tsv  | awk -vOFS="\t" '{print $1,$2,$3}' > ${INFOFILE}
SNUM=$( wc -l ${INFOFILE} | awk '{print $1}' )

echo "[STATUS] Set up infofile with $SNUM pairs using table $(basename ${MIXOBS_TAB})"

# --------------------------
# 1. Run command to extract:
# --------------------------
qsub -cwd -t 1-$SNUM -P compbio_lab -tc 2500 -l h_vmem=8G -l h_rt=0:45:00 -N extract_bw_pk_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_bw_on_peaks.sh -i $INFOFILE -o $RESULTDIR -e $EXTEND -m $MDPT"

# ------------------------------
# 2. Aggregate stats from above:
#    Plot stats.
# ------------------------------
PLOTPREF=${IMGDIR}/bwpk_intersect/all_bwpk_H3K27ac_in_DNase-seq
mkdir -p ${IMGDIR}/bwpk_intersect/

cat $RESULTDIR/BSS*stats.tsv | awk 'NR == 1{print $0}$1 != "metric"{print $0}' > $RESULTDIR/all_bwpk_statistics.tsv

qsub -cwd -P compbio_lab -tc 2500 -l h_vmem=8G -l h_rt=0:45:00 -N plot_bwpk_aggregate_${UQID} -hold_jid extract_bw_pk_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "conda activate mv_env >> /dev/null; R --slave -f ${BINDIR}/plot_aggregate_stats_bwpk.R --args $RESULTDIR/all_bwpk_statistics.tsv 'H3K27ac in DNase-seq peaks' $PLOTPREF $INFOFILE; conda deactivate"

# 3. Aggregate plots
