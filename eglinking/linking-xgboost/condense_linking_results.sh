#!/usr/bin/env bash
# ------------------------
# Collate linking results:
# ------------------------
start=`date +%s`
hostname -f
# Arguments:
export DISTANCE="0"
export MARK="0"
export ROADMAP="0"

if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [-t COMMAND] [OPTIONS] 
    -d     Distance at which to boost with distance pred. [default: $DISTANCE]
    -m     Whether to add mark information for prediction (0/1) [default: $MARK]
    -r     Whether to use logistic reg. roadmap method (0/1) [default: $ROADMAP]"
    exit 1
fi

while getopts d:m:r: o
do      case "$o" in
    d)      export DISTANCE="$OPTARG";;
    m)      export MARK="$OPTARG";;
    r)      export ROADMAP="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

# Load in variables:
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]; then
    source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
else
    source /home/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
fi
SRCDIR="$BINDIR/eglinking/linking-xgboost/linking/";
EGDIR="$BINDIR/eglinking/linking-xgboost/"
DATADIR="$DBDIR/linking_data/";

# Make directory:
PREDDIR=${DATADIR}predictions
suffix="pred"
if [[ "${ROADMAP}" == "1" ]]; then
    PREDDIR="${PREDDIR}_logistic"
    suffix="${suffix}_logistic"
else
    if [[ "${MARK}" == "1" ]]; then
        PREDDIR="${PREDDIR}_mark"
        suffix="${suffix}_mark"
    fi

    if [[ "$DISTANCE" != "0" ]]; then
        PREDDIR="${PREDDIR}_d${DISTANCE}"
        suffix="${suffix}_d${DISTANCE}"
    fi
fi
PREDDIR="${PREDDIR}/"

echo "[STATUS] Predictions are in ${PREDDIR}"
mkdir -p $PREDDIR $PREDDIR/collated
cd $PREDDIR

$TMP_DIR

MAPCOORD=${CORECOORD%%txt}.mapping.txt
awk -vOFS="\t" '{print $1"_"$2"_"$3, $4}' $CORECOORD | sort +0 -1 > $MAPCOORD

# ----------------------------
# Condense the collated files:
# ----------------------------
FULLFILE=${PREDDIR}/collated/allsamples_collated_${suffix}.tsv
while read -r sample; do
    echo $sample
    collfile=${PREDDIR}/collated/${sample}_collated_${suffix}.tsv.gz
    zcat $collfile | awk -vOFS="\t" '{print $1"_"$2"_"$3,$4,$5}' | sort +0 -1 | join -1 1 -2 1 - $MAPCOORD | awk -vOFS="\t" -v samp=$sample '{print $4,$2,$3,samp}' >> $FULLFILE
done < $DATADIR/Enhancer_matrix_names.txt

gzip $FULLFILE

