#!/usr/bin/env bash
# -------------------------------------
# Submit runs for enhancer-gene linking
# Basic XGBoost method:
# qsub -t 1-4998 -l h_rt=2:00:00 -l h_vmem=20G -cwd -N pred_eg_xg -j y -b y -V -r y -o $DBDIR/out/ "${EGDIR}/submit_linking.sh -d 0 -m 0"
# Distance boosted (50k and 100k):
# qsub -t 1-4998 -l h_rt=2:00:00 -l h_vmem=20G -cwd -N pred_eg_xg_50k -j y -b y -V -r y -o $DBDIR/out/ "${EGDIR}/submit_linking.sh -d 50000"
# qsub -t 1-4998 -l h_rt=2:00:00 -l h_vmem=20G -cwd -N pred_eg_xg_100k -j y -b y -V -r y -o $DBDIR/out/ "${EGDIR}/submit_linking.sh -d 100000"
# With the mark data:
# qsub -t 1-4998 -l h_rt=2:00:00 -l h_vmem=20G -cwd -N pred_eg_xg_mark -j y -b y -V -r y -o $DBDIR/out/ "${EGDIR}/submit_linking.sh -m 1"
# qsub -t 1-4998 -l h_rt=2:00:00 -l h_vmem=20G -cwd -N pred_eg_xg_mark -j y -b y -V -r y -o $DBDIR/out/ "${EGDIR}/submit_linking.sh -m 1 -d 50000"
# -------------------------------------
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
if [[ "${ROADMAP}" == "1" ]]; then
    PREDDIR="${PREDDIR}_logistic"
else
    if [[ "${MARK}" == "1" ]]; then
        PREDDIR="${PREDDIR}_mark"
    fi

    if [[ "$DISTANCE" != "0" ]]; then
        PREDDIR="${PREDDIR}_d${DISTANCE}"
    fi
fi
PREDDIR="${PREDDIR}/"

echo "[STATUS] Predictions will go to ${PREDDIR}"
mkdir -p $PREDDIR
cd $PREDDIR

# Make all possible predictions:
PREDTAB=${DATADIR}/predictions_table.tsv
if [[ ! -s $PREDTAB ]]; then
    join -j2 $DATADIR/Enhancer_matrix_names.txt $DATADIR/enh_chromhmm_states.txt > $PREDTAB
fi

SAMPLE=$( awk -v id=$SGE_TASK_ID 'NR==id{sub("^ ","",$0); print $1}' ${PREDTAB} )
CHROMHMM_STATE=$( awk -v id=$SGE_TASK_ID 'NR==id{sub("^ ","",$0); print $2}' ${PREDTAB} )
echo "[STATUS} Will predict for: ${SAMPLE} / ${CHROMHMM_STATE}"


cmd="$EGDIR/run_linking.sh -i $SAMPLE -s $CHROMHMM_STATE -d $DISTANCE -m $MARK -r $ROADMAP"
echo "$cmd"
bash -c "$cmd"


# ---------------------
# For selective re-run:
# ---------------------
if [[ "0" == "1" ]]; then
    for ITER in `seq 1 4998`; do 
        RUNPARAM=$( awk -v id=$ITER 'NR==id{sub("^ ","",$0); print $0}' ${PREDTAB} )
        BASENAME=`echo "${RUNPARAM}_2.5" | tr ' ' '_'`;
        if [[ ! -s ${BASENAME}.txt ]]; then
            echo "${ITER} ${BASENAME}"
            echo "${ITER} ${BASENAME}" >> to_rerun.tsv
            qsub -t ${ITER} -l h_rt=3:00:00 -l h_vmem=25G -cwd -N pred_eg_xg -j y -b y -V -r y -o $DBDIR/out/ "${EGDIR}/submit_linking.sh -d $DISTANCE -m $MARK -r $ROADMAP"
        else 
            echo "$ITER OK"
        fi
    done
    cat to_rerun.tsv
    wc to_rerun.tsv


    echo "[STATUS] There are $(ls $PREDDIR/BSS*.txt | wc -l) link sets, of which $( find $PREDDIR/BSS*.txt -not -empty -ls | wc -l ) are not empty"

    # Otherwise run all serially:
    conda activate pytorch_env

    for ITER in `seq 1 4998`; do 
        SAMPLE=$( awk -v id=$ITER 'NR==id{sub("^ ","",$0); print $1}' ${PREDTAB} )
        CHROMHMM_STATE=$( awk -v id=$ITER 'NR==id{sub("^ ","",$0); print $2}' ${PREDTAB} )
        BASENAME=`echo "${SAMPLE}_${CHROMHMM_STATE}_2.5" | tr ' ' '_'`;
        if [[ ! -s ${BASENAME}.txt ]]; then
            echo "$ITER"
            echo "${ITER} ${BASENAME}" >> to_rerun.tsv
            echo "[STATUS} Will predict for: ${SAMPLE} / ${CHROMHMM_STATE}"
            cmd="$EGDIR/run_linking.sh -i $SAMPLE -s $CHROMHMM_STATE -d $DISTANCE -m $MARK -r $ROADMAP"
            echo "$cmd"
            bash -c "$cmd"
        else 
            echo "$ITER OK"
        fi
    done
fi

