#!/bin/bash
# ------------------------------------------------
# Get bigwig and transfer to the ENCODE DCC server
# ------------------------------------------------
# Run as: 
# qsub -cwd -t 3-833 -l h_vmem=2G -l h_rt=00:30:00 -N dcc_expt_chmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/encode_submission/run_encode_submission_chromHMM.sh"
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

export DCC_LAB="manolis-kellis"
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

for SGE_TASK_ID in `seq 2 834`; do
    # Arguments:
    TASK=${SGE_TASK_ID}
    echo "[STATUS] Running TASK: $TASK"
    ANNTSV=$ANNDIR/encode_chmm_annotation_submission_metatable.tsv
    FILETSV=$ANNDIR/encode_chmm_file_submission_metatable.tsv
    TMP_DIR=${TMP}/encode_submission_chmm_${TASK}_${RANDOM}
    mkdir -p ${TMP_DIR}
    cd $TMP_DIR

    # Make $TASK-specific reduced metadata files:
    ATSV=${TMP_DIR}/annotation_metadata.tsv
    FTSV=${TMP_DIR}/file_metadata.tsv

    awk -v task=$TASK -F"\t" -vOFS="\t" 'NR==1{print $0}NR==task{print $0}' $ANNTSV > $ATSV

    # Get BSSID:
    bid=$( awk 'NR==2{sub(".*BSS","BSS",$0); sub(": .*","", $0); print $0}' $ATSV )
    nrun=$( grep "$bid" $DBDIR/chmm_enc.tsv -c )
    if [[ "$nrun" == "0" ]]; then
        echo "$bid"
        echo "[STATUS] Annotation metadata TSV:"
        cat $ATSV

        # Alias (pull all files):
        echo "\n[STATUS] File metadata TSV:"
        annalias=$( awk -F"\t" 'NR==2{print $1}' $ATSV ) 
        awk -v dataset=$annalias -F"\t" -vOFS="\t" 'NR==1{print $0}$1 == dataset{print $0}' $FILETSV > $FTSV
        cat $FTSV

        # -----------------------------
        # Submission of the annotation:
        # -----------------------------
        # Submit the experiment annotations:
        echo "[STATUS] Submitting ANNOTATION"
        cmd="python ${EU_PATH} -m $MODE -p annotation -i $ATSV"
        echo "$cmd"
        bash -c "$cmd"

        # --------------------------------
        # Submission of each of the files:
        # --------------------------------
        # Get dataset name and URL:
        cd $TMP_DIR
        for LNUM in `seq 2 5`; do 
            echo $LNUM
            FNUMTSV=${TMP_DIR}/file_metadata_${LNUM}.tsv
            awk -v lnum=$LNUM -F"\t" -vOFS="\t" 'NR==1{print $0}NR==lnum{print $0}' $FTSV > $FNUMTSV
            filename=$( awk -F"\t" -vOFS="\t" 'NR==1{a=0; for(i=1; i<=NF; i++){if ($i == "submitted_file_name"){a=i}}}NR>1{print $a}' $FNUMTSV ) 
            assembly=$( awk -F"\t" -vOFS="\t" 'NR==1{a=0; for(i=1; i<=NF; i++){if ($i == "assembly"){a=i}}}NR>1{print $a}' $FNUMTSV ) 
            filetype=$( awk -F"\t" -vOFS="\t" 'NR==1{a=0; for(i=1; i<=NF; i++){if ($i == "file_format"){a=i}}}NR>1{print $a}' $FNUMTSV ) 
            echo "Name:$filename\tGenome: $assembly\tType: $filetype"

            if [[ "$assembly" == "GRCh38" ]]; then
                genome="hg38"
            else 
                genome="hg19"
            fi

            # NOTE: For now, only submit hg19:
            # if [[ "$genome" == "hg19" ]]; then
            # Copy from the correct directory:
            CPUBDIR=/web/personal/cboix/epimap/ChromHMM
            FDIR=$CPUBDIR/observed_aux_18_${genome}/CALLS/
            if [[ "$filetype" == "bigBed" ]]; then
                FDIR=$FDIR/bigBed/
            fi

            # Add file_format_type: # ALREADY THERE.
            # awk -F"\t" -vOFS="\t" 'NR==1{print $0, "file_format_type"}NR==2{print $0,"bed9"}' $FNUMTSV > $TMP_DIR/tmp_filenum.tsv
            # mv $TMP_DIR/tmp_filenum.tsv $FNUMTSV

            cp $FDIR/$filename .

            # Submit the experiment files:
            echo "[STATUS] Submitting FILE METADATA + FILE"
            cmd="python $EU_PATH -m $MODE -p file -i $FNUMTSV"
            echo "$cmd"
            bash -c "$cmd"

            rm $filename
            # fi
        done
    fi
    rm -rf ${TMP_DIR}
done

end=`date +%s`
runtime=$((end-start))
echo "Finished ChromHMM DCC upload in $runtime seconds."
