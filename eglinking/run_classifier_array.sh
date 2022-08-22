# Grid Engine options
#$ -N run_classifier
#$ -cwd
#$ -P compbio_lab
#$ -l h_vmem=30G
#$ -l h_rt=16:00:00
# #$ -tc 250
# #$ -M cboix@mit.edu 
# #$ -m a 
#$ -j y
#$ -b y 
#$ -V 
#$ -r y 
#$ -o /broad/compbio/cboix/EPIMAP_ANALYSIS/db/linking_data/out/
#$ -e /broad/compbio/cboix/EPIMAP_ANALYSIS/db/linking_data/out/
#$ -t 1-833
start=`date +%s`
hostname -f

# ----------
# Arguments:
# ----------
INFOFILE=0
STATE_LIST=0
THRESHOLD=1
if [ $# -eq 0 ]; then
    echo >&2 "USAGE: source $0 [OPTIONS] 
    -i  Infofile - list of samples to run [defaults to 833 samples]
    -d  Saved dataset [required]
    -s  Enhancers - ChromHMM state list [optional]
    -t  Threshold [default: $THRESHOLD]"
    exit 1
fi

while getopts i:d:s:t: o
do      case "$o" in
    i)		INFOFILE="$OPTARG";;
    d)		DATA_PICKLE="$OPTARG";;
    s)		STATE_LIST="$OPTARG";;
    t)      THRESHOLD="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done


# -------------------
# Directories/config:
# -------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

DATA_DIR="$DBDIR/linking_data"
SAMPLE_DIR="$DATA_DIR/sample"
SRC_DIR="$BINDIR/eglinking/src"

if [[ "$INFOFILE" == "0" ]]; then
    INFOFILE=$DATA_DIR/mark_matrix_names.txt 
fi

if [[ "$STATE_LIST" == "0" ]]; then
    STATE_LIST=${SAMPLE_DIR}/state_list.txt 
fi

# Task-specific arguments:
samp=$(awk -v id=${SGE_TASK_ID} 'NR==id{print $0}' $INFOFILE )
OUTFILE=$DATA_DIR/links/${samp}_t${THRESHOLD}.tsv

# Run classifier:
echo "[STATUS] Running on ${samp} to create outfile: $OUTFILE"
conda activate linking

$SRC_DIR/classifier.py --saved $DATA_PICKLE --enh_state_list $STATE_LIST --threshold $THRESHOLD --which_first $samp --output $OUTFILE

EXIT_CLASSIF=$?
echo "[STATUS] Finished run with exitcode ${EXIT_CLASSIF}"

end=`date +%s`
runtime=$((end-start))
echo "[STATUS] Finished classifier run for ${samp} sucessfully in $runtime seconds."
