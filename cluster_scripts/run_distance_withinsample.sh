#!/bin/bash
# ================================
# Run distance between all obs/imp
# in a single sample across marks.
# ================================
start=`date +%s`
hostname -f
TASK=${SGE_TASK_ID}
TMP_DIR=${TMP}/sampdist_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR}

if [[ $# -lt 1 ]];then
    INFOFILE=${SAMPLEMARK_TAB}
else
    INFOFILE=$1
fi

if [[ $# -lt 2 ]];then
    OUTDIR=${MARKDIST_DIR}
else
    OUTDIR=$2
fi

SAMPLE=$( cut -f1 ${INFOFILE} | sort -u | sed "${TASK}q;d" - | awk -v FS="\t" '{print $1}' )

echo "[STATUS] Running distances within sample: $SAMPLE"

# ---------------------------------
# Make sample specific samplesheet:
# NOTE: sample and mark are swapped 
# to run dist within the mark.
# ---------------------------------
TMP_TAB=${TMP_DIR}/specific_table.tsv
grep -e "^$SAMPLE" $ALL_UQ_TAB | awk -vOFS="\t" '{split($1,a,"_"); print $2"_"a[2], a[1],$3}' > ${TMP_TAB}

echo "[STATUS] Using the following samplesheet table:"
cat ${TMP_TAB}

while read imark isample ifile; do
    echo "[STATUS] Compute global distance between datasets with ComputeGlobalDist"
    # Run distance:
    if [[ ! -s ${OUTDIR}/${imark}_${isample}.txt ]]; then
        echo "[STATUS] Computing for $imark in ${isample}. File is $ifile"
        # NOTE: Mark and sample are swapped in order to evaluate within sample not within mark.
        java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} ComputeGlobalDist -s ${imark} ${isample} ${IMPUTED_DIR} ${TMP_TAB} ${CHROMSIZES_noY} ${OUTDIR}
    else
        echo "[STATUS] Dist for $imark in ${isample} exists at ${OUTDIR}/${imark}_${isample}.txt"
    fi
done < <(grep "_obs" ${TMP_TAB} )

rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished $COMMAND sucessfully in $runtime seconds."
