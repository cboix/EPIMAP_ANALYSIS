#!/bin/bash
# -----------------------------------
# Lift over files from hg19 to hg38
# For chromHMM models
# -----------------------------------
NSTATES=15
MNAME=""
if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [OPTIONS] 
    -i     Info table [required]
    -n     Number of states model (15, 18, 25) [default: $NSTATES]
    -m     MODEL to use (overrides -n) (OPTIONAL)"
    exit 1
fi

while getopts n:m:i: o
do      case "$o" in
    i)      INFOFILE="$OPTARG";;
    n)		NSTATES="$OPTARG";;
    m)		MNAME="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

# 0. Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh ${NSTATES}
if [[ "$MNAME" != "" ]];then
    MODELPREF=$MNAME
else
    MODELPREF=$MODELNAME
fi

SAMPLE=$( cut -f1 $INFOFILE | sort -u | sed "${SGE_TASK_ID}q;d" )

echo "[STATUS] OPTIONS for liftOver set to:"
echo "[STATUS] nstates $NSTATES model $MODELPREF sample $SAMPLE"

export INCALLDIR=${CALLDIR}/${MODELPREF}
export OUTCALLDIR=${CALLDIR}/${MODELPREF}_hg38
mkdir -p ${OUTCALLDIR}/STATEBYLINE

INSEG=${INCALLDIR}/${SAMPLE}_${NSTATES}_CALLS_segments.bed.gz
OUTSEG=${OUTCALLDIR}/${SAMPLE}_${NSTATES}_CALLS_segments.bed 
TMPMAP=${OUTSEG%.bed}_unmapped

if [[ ! -s ${OUTSEG}.gz ]] || [[ ${INSEG} -nt ${OUTSEG}.gz ]]; then
    echo "[STATUS] Lifting over $INSEG"
    liftOver ${INSEG} ${LOCHAIN} ${OUTSEG} ${TMPMAP}
    gzip -f $OUTSEG

    rm ${TMPMAP}
fi

# while read chr size; do
#     suffix=CALLS_PER_LINE_${chr}_statebyline.txt
#     INFILE=${INCALLDIR}/STATEBYLINE/${SAMPLE}_${NSTATES}_${suffix}.gz
#     OUTFILE=${OUTCALLDIR}/STATEBYLINE/${SAMPLE}_${NSTATES}_${suffix}
#     if [[ ! -s ${OUTFILE}.gz ]] || [[ ${INFILE} -nt ${OUTFILE}.gz ]]; then
#         echo $chr
#         liftOver $INFILE ${LOCHAIN} $OUTFILE ${OUTFILE}_unmapped
#         gzip -f $OUTFILE
#         rm ${OUTFILE}_unmapped
#     fi
# done < ${CHROMSIZES_noY}
