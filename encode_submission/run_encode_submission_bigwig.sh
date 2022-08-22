#!/bin/bash
# ------------------------------------------------
# Get bigwig and transfer to the ENCODE DCC server
# Updated 05/11/21
# ------------------------------------------------
# Run as: 
# qsub -cwd -t 11-14909 -tc 500 -l h_vmem=2G -l h_rt=00:30:00 -N dcc_expt_bw -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/encode_submission/run_encode_submission_bigwig.sh"
# qsub -cwd -t 2-14909 -tc 500 -l h_vmem=3G -l h_rt=01:00:00 -N dcc_expt_bw -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/encode_submission/run_encode_submission_bigwig.sh"
# Directories/Config:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh
start=`date +%s`
hostname -f

# Required for submission:
# yes | conda install -c conda-forge requests awscli google-api-python-client inflection jsonschema urllib3 -n pytorch_env 
conda activate pytorch_env

# NOTE: Must have these somewhere:
# export DCC_API_KEY='xxxxxx'
# export DCC_SECRET_KEY='xxxxxxxx'

export DCC_API_KEY="R7GOMSZL"
export DCC_SECRET_KEY="4ftc3mvvllhkdpe2"

# Using eu_register:
export EU_RELEASE=$SFTDIR/encode_utils
# Update your PYTHONPATH environment variable as follows:
export PYTHONPATH=${EU_RELEASE}:${PYTHONPATH}
# Update your PATH environment variable as follows:
export PATH=${EU_RELEASE}/encode_utils/MetaDataRegistration:${PATH}
export PATH=${EU_RELEASE}/encode_utils/scripts:${PATH}
export EU_PATH=${EU_RELEASE}/encode_utils/MetaDataRegistration/eu_register.py

# For running submission:
MODE=prod
# MODE=dev

# for SGE_TASK_ID in `seq 2 14909`; do
# Arguments:
TASK=${SGE_TASK_ID}
echo "[STATUS] Running TASK: $TASK"
ANNTSV=$ANNDIR/encode_bigwig_annotation_submission_metatable.tsv
FILETSV=$ANNDIR/encode_bigwig_file_submission_metatable.tsv
TMP_DIR=${TMP}/encode_submission_bw_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR}

# Make $TASK-specific reduced metadata files:
ATSV=${TMP_DIR}/annotation_metadata.tsv
FTSV=${TMP_DIR}/file_metadata.tsv

# Pull TSVs:
awk -v task=$TASK -F"\t" -vOFS="\t" 'NR==1{print $0}NR==task{print $0}' $ANNTSV > $ATSV
awk -v task=$TASK -F"\t" -vOFS="\t" 'NR==1{print $0}NR==task{print $0}' $FILETSV > $FTSV

# Check file transfer completion:
fid=$( awk -F"\t" 'NR==2{print $7}' $FTSV )
nrun=$( grep "$fid" $DBDIR/bw_enc.tsv -c )

if [[ "$nrun" == "0" ]]; then
    # Prune targets from file:
    assay=$( awk -F"\t" -vOFS="\t" 'NR==1{a=0; for(i=1; i<=NF; i++){if ($i == "assay_term_name"){a=i}}}NR==2{print $a}' $ATSV ) 
    if [[ "$assay" == "ATAC-seq" ]] || [[ "$assay" == "DNase-seq" ]]; then
        awk -F"\t" -vOFS="\t" 'NR==1{a=0; printf $1; for(i=2; i<=NF; i++){if ($i == "targets"){a=i} else {printf "\t"$i }}; printf "\n"}NR==2{ printf $1; for(i=2; i<=NF; i++){if (i != a){ printf "\t"$i }}; printf "\n" }' $ATSV > tmp
        mv tmp $ATSV
    fi
    echo "[STATUS] Annotation metadata TSV:"
    cat $ATSV

    echo "\n[STATUS] File metadata TSV:"
    cat $FTSV

    # Get dataset name and URL:
    filename=$( awk -v task=$TASK -F"\t" -vOFS="\t" 'NR==1{a=0; for(i=1; i<=NF; i++){if ($i == "submitted_file_name"){a=i}}}NR==2{print $a}' $FTSV ) 
    datatype=$( awk -v task=$TASK -F"\t" -vOFS="\t" 'NR==2{if ($1 ~ /imputed/){a = "imputed"} else {a = "observed"}; print a}' $FTSV ) 
    fileurl="https://epigenome.wustl.edu/epimap/data/${datatype}/${filename}"
    echo "Name: $filename\tSet: $datatype"

    # Pull the data into temporary directory:
    echo "[STATUS] Pulling data from ${fileurl}"
    cd $TMP_DIR
    wget $fileurl --quiet
    echo "$( ls -shclt ${TMP_DIR} )"

    # Test successful pull:
    # gzip -t $filename
    checkcmd="import pyBigWig; bw = pyBigWig.open('"$filename"'); print(bw.isBigWig());" 
    checkvalid=$( python -c "$checkcmd" )

    if [[ "$checkvalid" == "True" ]]; then
        # Submit the experiment annotations:
        echo "[STATUS] Submitting ANNOTATION"
        cmd="python ${EU_PATH} -m $MODE -p annotation -i $ATSV"
        echo "$cmd"
        bash -c "$cmd"

        # Submit the experiment files:
        echo "[STATUS] Submitting FILE METADATA + FILE"
        cmd="python ${EU_PATH} -m $MODE -p file -i $FTSV"
        echo "$cmd"
        bash -c "$cmd"
    else 
        echo "[STATUS] Failed getting file, validity check result is $checkcmd"
    fi

    rm $filename
    rm -rf ${TMP_DIR}
fi

end=`date +%s`
runtime=$((end-start))
echo "Finished bigWig DCC upload in $runtime seconds."
