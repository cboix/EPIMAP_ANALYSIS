#!/bin/bash
# -----------------------------------
# Output the signal distr. statistics
# For a specific track
# -----------------------------------
start=`date +%s`
hostname -f

if [[ $# -lt 1 ]];then
    INFOFILE=${ALL_UQ_TAB}
else
    INFOFILE=$1
fi

if [[ $# -lt 2 ]];then
    OUTDIR=${DISTR_DIR}
else
    OUTDIR=$2
fi

if [[ $# -lt 3 ]];then
    TASK=${SGE_TASK_ID}
else
    TASK=$3
fi

# Directories:
TMP_DIR=${TMP}/signaldistr_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR}

SAMPLE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
MARK=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $2}' )
FPREF=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $3}' )
echo "[STATUS] Getting distribution of $SAMPLE in $MARK"

# Output distr with header as prefix
if [[ "$( echo "$FPREF" | grep -e "^impute" - -c )" == "1" ]]; then
    DATADIR=$IMPUTED_DIR
else
    DATADIR=$CONVERTED_DIR
fi

# All counts, agnostic of state location
MAINTABLE=${OUTDIR}/${FPREF}_signalcounts.tsv
if [[ ! -s $MAINTABLE ]]; then
    # Read all files for relevant chr:
    while read chr size; do
        FILE=${DATADIR}/${chr}_${FPREF}.wig.gz 
        zcat $FILE | awk -vOFS="\t" 'NR==1{split("",a,"")}NR > 2{a[$1]++}END{for (i in a) print i,a[i]}' | sort -n >> ${OUTDIR}/${FPREF}_tmp.tsv
    done < ${CHROMSIZES_noY}

    # Process all together:
    awk -vOFS="\t" 'NR==1{split("",a,"")}NR > 2{a[$1] = a[$1] + $2}END{for (i in a) print i,a[i]}' ${OUTDIR}/${FPREF}_tmp.tsv | sort -n > $MAINTABLE 

    rm ${OUTDIR}/${FPREF}_tmp.tsv
fi

STATETABLE=${OUTDIR}/${FPREF}_signal_bystate.tsv
if [[ ! -s ${OUTFILE} ]]; then
    source activate pytorch_env
    python $BINDIR/extract_signaldistr_statespecific.py main --maintrack ${FPREF} --output ${STATETABLE} --mark $MARK --trackdir ${DATADIR}/
    source deactivate
fi

rm -rf $TMP_DIR

end=`date +%s`
runtime=$((end-start))
echo "Finished getting distribution in $runtime seconds."
