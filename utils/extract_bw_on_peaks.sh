#!/bin/bash
# -----------------------------------
# Extract BW around peak midpoints
# -----------------------------------
start=`date +%s`
hostname -f
# Default arguments:
TASK=0
MDPT=0
EXTEND=250

if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [OPTIONS] 
    -i     Infofile (id peakfile bwfile) [required]
    -o     Output directory [required]
    -e     +/- bp extension from midpt (OPTIONAL) [default: $EXTEND]
    -m     If not 0, midpt column in peak file (OPTIONAL) [default: $MDPT]
    -t     Task number (OPTIONAL) [default: SGE_TASK_ID]"
    exit 1
fi

while getopts i:o:m:e:t: o
do      case "$o" in
    i)		INFOFILE="$OPTARG";;
    o)		RESULTDIR="$OPTARG";;
    e)		EXTEND="$OPTARG";;
    m)		MDPT="$OPTARG";;
    t)		TASK="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

if [[ "$TASK" == "0" ]]; then
    TASK=${SGE_TASK_ID}
fi

SAMPLE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
PEAKSFILE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $2}' )
BWFILE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $3}' )
FINALFILE=${RESULTDIR}/${SAMPLE}_pm${EXTEND}.tsv.gz
FINALPLOT=${RESULTDIR}/${SAMPLE}_pm${EXTEND}_heatmap.png
FINALSTATS=${RESULTDIR}/${SAMPLE}_pm${EXTEND}_stats.tsv
echo "[STATUS] Extracting peak regions $PEAKSFILE from bw $BWFILE for sample id $SAMPLE"

TMP_DIR=${TMP}/extract_bw_peaks_${SAMPLE}_${RANDOM}
mkdir -p ${TMP_DIR}

# Identify where files come from:
IMPPK=$( echo "$PEAKSFILE" | grep "impute" | wc -l )
if [[ "$IMPPK" == "1" ]];then 
    PKDIR=${PKIMP_DIR} 
else
    PKDIR=${PKOBS_DIR}
fi

IMPBW=$( echo "$BWFILE" | grep "impute" | wc -l )
if [[ "$IMPBW" == "1" ]];then 
    BWDIR=${IMPUTED_DIR} 
else
    BWDIR=${CONVERTED_DIR}
fi

# ---------------------
# Perform intersection:
# ---------------------
if [[ ! -s $FINALFILE ]] || [[ "$( zcat $FINALFILE | wc -l )" == "0" ]]; then
    echo "[STATUS] Performing intersection"
    # Make temporary bedfile:
    BEDFILE=${TMP_DIR}/regions.bed
    if [[ "$MDPT" == "0" ]]; then 
        echo "[STATUS] Creating bedfile by taking midpoint of start + end"
        zcat $PKDIR/$PEAKSFILE | awk -F"\t" -vOFS="\t" -v md=$MDPT -v ext=$EXTEND '{md=($2+$3)/2; s=(md-ext>0?md-ext:0); print $1, s, md+ext, NR}' > $BEDFILE
    else
        echo "[STATUS] Creating bedfile with specified midpoint (column #$MDPT)"
        zcat $PKDIR/$PEAKSFILE | awk -F"\t" -vOFS="\t" -v md=$MDPT -v ext=$EXTEND '{s=($md-ext>0?$md-ext:0); print $1, s, $md+ext, NR}' > $BEDFILE
    fi

    NREGIONS=$( cat $BEDFILE | wc -l )
    echo "[STATUS] There are $NREGIONS regions for intersection"

    # NOTE: For BW, have multiple chromosomes.
    TMPINT=${TMP_DIR}/intersection.tsv
    rm $TMPINT
    while read chr size; do 
        echo "[STATUS] Processing $chr"
        BWCHROM=${BWDIR}/${chr}_${BWFILE}.wig.gz
        # Turn into bedfile:
        echo "[STATUS] Making bedfile:"
        zcat $BWCHROM | awk -vOFS="\t" -v chr=$chr 'NR > 2{print chr, 1 + 25*(NR-3), 25*(NR-2),$1}' | gzip -c > ${TMP_DIR}/tmp_bw.bed.gz
        # Intersect temporary file with bedfile:
        echo "[STATUS] Intersecting:"
        bedtools intersect -a $BEDFILE -b ${TMP_DIR}/tmp_bw.bed.gz -wa -wb | awk -vOFS="\t" '{m=($2+$3)/2; s=($2>$6?$2:$6); e=($3<$7?$3:$7); print $1"_"m, s-m,e-m,$8}' >> ${TMPINT}
    done < ${CHROMSIZES_noY}

    gzip -c $TMPINT > $FINALFILE

    echo "[STATUS] Intersection is in ${FINALFILE}"
    ls -sh $FINALFILE
else
    echo "[STATUS] File exists and can be found at ${FINALFILE}"
    ls -sh $FINALFILE
fi

# --------------------------------
# Intersection statistics + plots:
# --------------------------------
if [[ ! -s $FINALSTATS ]] || [[ ! -s $FINALPLOT ]]; then
    echo "[STATUS] Getting intersection statistics and plotting intersection"
    conda deactivate >> /dev/null  # In case of loaded env.
    source activate mv_env >> /dev/null
    PLOTCMD="R --slave -f ${BINDIR}/plot_stats_bw_pk_intersection.R --args $FINALFILE $SAMPLE $FINALSTATS $FINALPLOT"
    echo "Running:\n${PLOTCMD}"
    bash -c "${PLOTCMD}"
    conda deactivate >> /dev/null
fi

rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished intersection sucessfully in $runtime seconds."
