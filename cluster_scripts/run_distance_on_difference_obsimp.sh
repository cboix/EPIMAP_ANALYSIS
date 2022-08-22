#!/bin/bash
# --------------------------------------------
# Tool to run/parallelize difference analyses:
# 
# Difference calculations:
# DIFFERENCE = Raw difference of two tracks
# SCALEDIFF = Rescaled diff. of two tracks
# 
# For distance, 
# DISTANCE = Within-sample, cross-mark
# SAMPDIST = Cross-sample, within-mark
# --------------------------------------------
if [[ $# -lt 1 ]]; then
    echo "USAGE: $(basename $0) [COMMAND] [INFOFILE] (optional [TASK])" >&2
    echo '  [COMMAND]: Command to execute out of:' >&2
    echo '             DIFFERENCE, SCALEDIFF, DISTANCE, SAMPDIST' >&2
    exit 1
fi
COMMAND=$1

if [[ $# -gt 1 ]]; then 
    INFOFILE=$2
else
    INFOFILE=$SAMPLEMARK_TAB
fi

start=`date +%s`
hostname -f

TASK=${SGE_TASK_ID}
TMP_DIR=${TMP}/diffdist_${COMMAND}_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR}

# For distance:
TMP_TAB=${TMP_DIR}/specific_table.tsv

# Get sample/mark from table (available or to_impute):
SAMPLE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
MARK=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $2}' )
OBS_FILE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $3}' )
IMP_FILE="impute_${SAMPLE}_${MARK}"
DIFF_FILE="difference_${SAMPLE}_${MARK}"

echo "[STATUS] Running command $COMMAND for observed vs imputed $SAMPLE and $MARK"

if [[ "$COMMAND" == "DIFFERENCE" ]]; then
    # Raw difference:
    TMPFILE=${TMP_DIR}/tmpfile.tsv
    while read chr size; do
        # Collate, normalize the observed side by chromosome, and subtract (cap at 0)
        zcat ${CONVERTED_DIR}/${chr}_${OBS_FILE}.wig.gz > ${TMPFILE}
        zcat ${IMPUTED_DIR}/${chr}_${IMP_FILE}.wig.gz | pr -mts ${TMPFILE} - | awk -v name="${SAMPLE}_${MARK}_difference" -v chr=$chr 'NR==1{print "track type=wiggle_0 name="name}NR==2{print "fixedStep chrom="chr" start=1 step=25 span=25"}NR >2{a=$1-$2; a=(a<0?0:a); print a}' | gzip -c > ${IMPUTED_DIR}/${chr}_${DIFF_FILE}.wig.gz
    done < ${CHROMSIZES_noY}

elif [[ "$COMMAND" == "SCALEDIFF" ]]; then
    # Regression to rescale, followed by difference:
    source activate pytorch_env
    python $BINDIR/calculate_diff_of_datasets.py main --observed ${OBS_FILE} --imputed ${IMP_FILE} --output ${DIFF_FILE}
    source deactivate

elif [[ "$COMMAND" == "DISTANCE" ]]; then
    # Make within-sample samplesheet:
    # NOTE: sample and mark are swapped to run dist within the mark.
    grep -e "^$SAMPLE" $ALL_TRACKS_TAB | awk -vOFS="\t" '$3 ~ /impute/{print $2, $1,$3}' > ${TMP_TAB}
    echo -e "diff_${MARK}\\t${SAMPLE}\\t${DIFF_FILE}" >> ${TMP_TAB}
    echo "[STATUS] Calculating cross-mark, within-sample distance. Using the following samplesheet table:"
    cat ${TMP_TAB}
    if [[ ! -s ${DIFFDIST_DIR}/diff_${MARK}_${SAMPLE}.txt ]]; then
        echo "[STATUS] Compute global distance between datasets with ComputeGlobalDist"
        # NOTE: Mark and sample are swapped in order to evaluate within sample not within mark.
        java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} ComputeGlobalDist -s "diff_${MARK}" ${SAMPLE} ${IMPUTED_DIR} ${TMP_TAB} ${CHROMSIZES_noY} ${DIFFDIST_DIR}
    fi

elif [[ "$COMMAND" == "SAMPDIST" ]]; then
    # Make within-mark samplesheet:
    grep "$MARK" $ALL_TRACKS_TAB | awk -vOFS="\t" '$3 ~ /impute/{print $1, $2, $3}' > ${TMP_TAB}
    echo -e "sampdiff_${SAMPLE}\\t${MARK}\\t${DIFF_FILE}" >> ${TMP_TAB}
    echo "[STATUS] Calculating within-mark, cross-sample distance. Using the following samplesheet table:"
    cat ${TMP_TAB}

    # Make samplesheet:
    if [[ ! -s ${DIFFDIST_DIR}/sampdiff_${SAMPLE}_${MARK}.txt ]]; then
        echo "[STATUS] Compute global distance between datasets with ComputeGlobalDist"
        # NOTE: Mark and sample are swapped in order to evaluate within sample not within mark.
        java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} ComputeGlobalDist -s "sampdiff_${SAMPLE}" ${MARK} ${IMPUTED_DIR} ${TMP_TAB} ${CHROMSIZES_noY} ${DIFFDIST_DIR}
    fi
else
    echo "Unknown impute difference command: $COMMAND"
    echo "Please use one of: DIFFERENCE, SCALEDIFF, DISTANCE, SAMPDIST"
    exit
fi

rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished $COMMAND sucessfully in $runtime seconds."
