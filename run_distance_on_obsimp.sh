#!/bin/bash
# ---------------------------------
# Run distance between all obs/imp:
# Modes:
# GWDIST (ChromImpute)
# Metrics within
# ---------------------------------
start=`date +%s`
hostname -f
if [[ $# -lt 1 ]]; then
    echo "USAGE: $(basename $0) [COMMAND] [INFOFILE] (optional [OUTDIR])" >&2
    echo '  [COMMAND]: Command to execute out of:' >&2
    echo '             GWDIST, PREPROCESS,' >&2
    echo '             FIXEDDIST, ALLDIST' >&2
    exit 1
fi
COMMAND=$1

if [[ $# -lt 2 ]];then
    INFOFILE=${ALL_UQ_TAB}
else
    INFOFILE=$2
fi

if [[ $# -lt 3 ]];then
    OUTDIR=${IMPDIST_DIR}
else
    OUTDIR=$3
fi

# NOTE: Set NREGIONS to negative vals to deactivate
if [[ $# -lt 4 ]];then
    NREGIONS=20000
else
    NREGIONS=$4
fi

TASK=${SGE_TASK_ID}
TMP_DIR=${TMP}/alldist_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR}

# Get sample/mark from table (available or to_impute):
SAMPLE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
MARK=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $2}' )
INPUTPREF=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $3}' )
echo "[STATUS] ChromImpute: Running full distance on imputed $SAMPLE and $MARK"

# Compare to all tracks in same mark:
TRACKSFILE=${TMP_DIR}/tracks_comparison.tsv
awk -v mark=$MARK '$2 == mark{print $3}' $ALL_UQ_TAB > $TRACKSFILE

if [[ "$COMMAND" == "GWDIST" ]]; then
    if [[ ! -s ${OUTDIR}/${SAMPLE}_${MARK}.txt ]]; then
        echo "[STATUS] Compute global distance between datasets with ComputeGlobalDist"
        java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} ComputeGlobalDist -s ${SAMPLE} ${MARK} ${IMPUTED_DIR} ${INFOFILE} ${CHROMSIZES} ${OUTDIR}
    fi

elif [[ "$COMMAND" == "PREPROCESS" ]]; then
    echo "[STATUS] Preprocess $INPUTPREF into subsetted cp file."
    source activate pytorch_env; 
    python ${BINDIR}/calculate_distance_of_datasets.py preprocess --maintrack ${INPUTPREF} --tracksfile ${TRACKSFILE} --output $OUTFILE --mark $MARK --fixed_region --nregions $NREGIONS
    source deactivate 

elif [[ "$COMMAND" == "FIXEDDIST" ]]; then
    echo "[STATUS] Compute distance between datasets in fixed regions (top $NREGIONS)"
    OUTFILE=${OUTDIR}/fixeddist_${NREGIONS}_${SAMPLE}_${MARK}.tsv
    if [[ ! -s $OUTFILE ]]; then
        source activate pytorch_env; 
        python ${BINDIR}/calculate_distance_of_datasets.py main --maintrack ${INPUTPREF} --tracksfile ${TRACKSFILE} --output $OUTFILE --mark $MARK --fixed_region --nregions $NREGIONS
        source deactivate 
    fi

elif [[ "$COMMAND" == "ALLDIST" ]]; then
    OUTFILE=${OUTDIR}/alldist_${NREGIONS}_${SAMPLE}_${MARK}.tsv
    if [[ ! -s $OUTFILE ]]; then
        source activate pytorch_env; 
        python ${BINDIR}/calculate_distance_of_datasets.py main --maintrack ${INPUTPREF} --tracksfile ${TRACKSFILE} --output $OUTFILE --mark $MARK --nofixed_region --nregions $NREGIONS
        source deactivate 
    fi

else
    echo "[STATUS] $COMMAND is not a valid command"
fi


rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished $COMMAND sucessfully in $runtime seconds."
