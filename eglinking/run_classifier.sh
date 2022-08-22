#!/usr/bin/env bash
# --------------------------
# Run classifier for linking
# Author: Benjamin T. James
# --------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

# Arguments:
SRC_DIR="$BINDIR/eglinking/src"


[[ $- == *i* ]] || die "Not an interactive shell, rerun with bash -i"
test -d $SRC_DIR || die "Source dir $SRC_DIR must exist";

# Create the appropriate conda environment:
export PATH="$SRC_DIR:$PATH";
conda activate linking

# Directories for linking:
DATA_DIR="$DBDIR/linking_data";
SAMPLE_DIR="$DATA_DIR/sample";

mkdir -p $DATA_DIR $DATA_DIR/out $DATA_DIR/links/
cd $DATA_DIR

# Arguments:
THRESHOLD=1.5
SAMPLE_LIST=${SAMPLE_DIR}/sample_list.txt 
STATE_LIST=${SAMPLE_DIR}/state_list.txt 

# Create statelist and sample list:
if [[ ! -s $SAMPLE_LIST ]]; then
    awk '{print $1}' ${SAMPLE_DIR}/shared_samples.txt > $SAMPLE_LIST
fi

if [[ ! -s $STATE_LIST ]]; then
    stlist=(E7 E8 E9 E10 E11 E15)
    rm $STATE_LIST
    for st in ${stlist[@]}; do
        echo $st >> $STATE_LIST
    done
fi

# DATA_PICKLE="$DATA_DIR/2020_03_06_18:43:17.pickle.gz"
DATA_PICKLE="$DATA_DIR/2020_03_31_22:06:21.pickle.gz"
if [[ ! -s $DATA_PICKLE ]]; then 
    # Run once for all 
    # $SRC_DIR/classifier.py $DATA_DIR $SAMPLE_DIR $SAMPLE_LIST $STATE_LIST $THRESHOLD
    $SRC_DIR/classifier.py --data_dir $DATA_DIR --sample_dir $SAMPLE_DIR --sample_list $SAMPLE_LIST --enh_state_list $STATE_LIST --threshold $THRESHOLD
    # NOTE: DATA PICKLE NAME WILL JUST BE THE TIME STAMP
fi

# Testing/debugging::
if [[ "1" == "0" ]]; then
    samp="BSS00004"
    # Old run style:
    # $SRC_DIR/classifier.py $DATA_PICKLE $STATE_LIST $THRESHOLD $samp
    OUTFILE=$DATA_DIR/links/${samp}_t${THRESHOLD}_tmp.tsv
    $SRC_DIR/classifier.py --saved $DATA_PICKLE --enh_state_list $STATE_LIST --threshold $THRESHOLD --which_first $samp --output $OUTFILE
fi

# Run as array across 833 samples:
INFOFILE=$DATA_DIR/mark_matrix_names.txt
BNUM=$(wc -l $INFOFILE | awk '{print $1}')
qsub -t 1-200 -l h_rt=16:00:00 $BINDIR/eglinking/run_classifier_array.sh -i $INFOFILE -d $DATA_PICKLE -s $STATE_LIST -t $THRESHOLD
# qsub -t 1-$BNUM $BINDIR/eglinking/run_classifier_array.sh -i $INFOFILE -d $DATA_PICKLE -s $STATE_LIST -t $THRESHOLD




