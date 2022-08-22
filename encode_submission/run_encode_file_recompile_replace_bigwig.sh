#!/bin/bash
# ------------------------------------------------
# Get bigwig and transfer to the ENCODE DCC server
# Updated 05/11/21
# ------------------------------------------------
# Run as: 
# qsub -cwd -t 1-223 -tc 500 -l h_vmem=3G -l h_rt=00:45:00 -N dcc_expt_recomp_bw -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/encode_submission/run_encode_file_recompile_replace_bigwig.sh"
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
export HLSCRIPT="$BINDIR/encode_submission/encode_file_replacement.py"

# For running submission:
MODE=prod
# MODE=dev

# for SGE_TASK_ID in `seq 1 223`; do 
# Arguments:
TASK=${SGE_TASK_ID}
echo "[STATUS] Running TASK: $TASK"
ANNTSV=$ANNDIR/encode_bigwig_annotation_submission_metatable.tsv
FILETSV=$ANNDIR/encode_bigwig_file_submission_metatable.tsv
TMP_DIR=${TMP}/encode_submission_bw_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR}

# Get the 
FID=$( awk -v task=$TASK -F"\t" 'NR==task{print $1}' $DBDIR/bw_replace.tsv )
FILENAME=$( awk -v task=$TASK -F"\t" 'NR==task{print $2}' $DBDIR/bw_replace.tsv )

# Check if file is in dblist:
nrun=$( grep "$FILENAME" $FILETSV -c )
# nrun=$( grep "$FILENAME" $FILETSV -n )

if [[ "$nrun" == "1" ]]; then
    cd $TMP_DIR
    # Create the bigWig file:
    INFOFILE=$PUBLIC_METADIR/all_released_tracks.tsv
    nline=$( grep "${FILENAME%%.bigWig}" $INFOFILE -n )
    num=${nline%%:*}

    # Get sample/mark from table (available or to_impute):
    SAMPLE=$( sed "${num}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
    MARK=$( sed "${num}q;d" ${INFOFILE} | awk -v FS="\t" '{print $2}' )
    FILE=$( sed "${num}q;d" ${INFOFILE} | awk -v FS="\t" '{print $3}' )
    echo "[STATUS] BDG to BW for $SAMPLE $MARK $FILE"

    # NOTE: Put each into its own mark directory:
    if [[ "$( echo ${FILE} | grep "impute" - -c )" == "1" ]]; then
        bdgdir=${IMPUTED_DIR}
        bdgclass="imputed"
    else
        bdgdir=${CONVERTED_DIR}
        bdgclass="observed"
    fi
    BWDIR=${PUBLIC_BWDIR}/${bdgclass}/${MARK}
    BDGDIR=${PUBLIC_BDGDIR}/${bdgclass}/${MARK}
    mkdir -p $BWDIR $BDGDIR

    tmpwig=${TMP_DIR}/${FILE}.wig
    tmpbw=${TMP_DIR}/${FILE}.bigWig
    rm $tmpwig
    while read chr size; do 
        echo ${chr}
        fname=${chr}_${FILE}.wig.gz
        cfile=${bdgdir}/$fname
        zcat $cfile | awk 'NR > 1' >> $tmpwig
        # Copy the BW file to the public bedgraph directory as well.
    done < ${CHROMSIZES_noY}

    cmd="wigToBigWig -clip $tmpwig ${CHROMSIZES_noY} $tmpbw"
    echo "$cmd"
    bash -c "$cmd"
    echo "$( ls -shclt ${TMP_DIR} )"

    # ----------------------------------
    # Replace the file at ENCODE portal:
    # ----------------------------------
    # Test successful creation:
    checkcmd="import pyBigWig; bw = pyBigWig.open('"$FILENAME"'); print(bw.isBigWig());" 
    checkvalid=$( python -c "$checkcmd" )
    echo "[STATUS] Validity check outputs: $checkvalid"

    if [[ "$checkvalid" == "True" ]]; then
        # Replace the file:
        # echo "[STATUS] Replacing FILE"
        # cmd="python $HLSCRIPT upload --fileid $FID --filepath $FILENAME --mode $MODE"
        echo "[STATUS] Checksum:"
        cmd="md5sum $FILENAME"
        echo "$cmd"
        bash -c "$cmd"
    else 
        echo "[STATUS] Failed getting file, validity check result is $checkcmd"
    fi

    rm $FILENAME
    rm -rf ${TMP_DIR}
fi

# done

end=`date +%s`
runtime=$((end-start))
echo "Finished bigWig DCC file replacement in $runtime seconds."
