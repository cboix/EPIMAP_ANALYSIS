#!/bin/bash
# ------------------------
# GO Enrichment submission
# TODO: For regions / BKG above the limit, subsample down
# ------------------------
# 0. Load in variables:
BKGFILE="NONE"
TASK=${SGE_TASK_ID}
if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $(basename $0) -i [INFOFILE] -o [RESULTDIR] (OTHERARGS)
    -i  Infofile [REQUIRED] 
            Must include BED Filename and UID for the dataset
    -o  Output directory for results [REQUIRED]
    -b  Background file (OPTIONAL)
    -t  Task (OPTIONAL)
            line from file to run. Otherwise runs SGE_TASK_ID: $SGE_TASK_ID"
    exit 1
fi

while getopts i:o:b:t: o
do      case "$o" in
    i)		INFOFILE="$OPTARG";;
    o)      RESULTDIR="$OPTARG";;
    b)		BKGFILE="$OPTARG";;
    t)		TASK="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

echo "[STATUS] INFOFILE: $INFOFILE"
echo "[STATUS] RESULTDIR: $RESULTDIR"
echo "[STATUS] BACKGROUND: $BKGFILE"
echo "[STATUS] TASK: $TASK"

# NOTE: Load after case catches variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
start=`date +%s` # For timing the pipeline

# Directories 
PUBROOT=/web/personal/${USER}
TMP_DIR=run_GREAT_analysis_${RANDOM}
mkdir -p ${PUBROOT}/${TMP_DIR} ${RESULTDIR}

# File attributes:
file=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
file_base=`basename ${file}`
TMP_FILE="${TMP_DIR}/regions_${RANDOM}.bed"
TMP_BKG="${TMP_DIR}/background_regions.bed"
file_new=${file_base%%.gz}_results

# NOTE: Label each peak by peak NR and tasknum:
# (Otherwise, get huge results)
if [[ $file =~ \.gz$ ]]; then
    zcat ${file} | awk -vOFS="\t" -v task=$TASK '$1 ~ /chr[0-9X]/{print $1, $2, $3, task"_"NR}' > ${PUBROOT}/${TMP_FILE}
else
    awk -vOFS="\t" -v task=$TASK '$1 ~ /chr[0-9X]/{print $1, $2, $3, task"_"NR}' ${file} > ${PUBROOT}/${TMP_FILE}
fi

# Add background if given:
if [[ "$BKGFILE" != "NONE" ]]; then
    cp ${BKGFILE} ${PUBROOT}/${TMP_BKG}
fi

# Check file:
HEADER=$( zcat ${RESULTDIR}/${file_new}.gz | awk -F"\t" 'NR==2{print $1}')
if [[ "$HEADER" != "# GREAT version 3.0.0" ]]; then
    rm ${RESULTDIR}/${file_new}.gz 
else
    echo "[STATUS] File exists and has correct format"
fi

# GREAT submission:
if [[ ! -s ${RESULTDIR}/${file_new}.gz ]] ||
    [[ ${file} -nt "${RESULTDIR}/${file_new}.gz" ]]; then
    echo "[STATUS] Processing $file with $( cat ${PUBROOT}/${TMP_FILE} | wc -l ) regions"

    # Properties of request:
    REQ_BASE="http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php"
    REQ_OUT="outputType=batch&requestSpecies=hg19"
    REQ_FILE="requestURL=http://www.broadinstitute.org/~${USER}/${TMP_FILE}"
    if [[ "$BKGFILE" != "NONE" ]]; then
        REQ_BKG="bgURL=http://www.broadinstitute.org/~${USER}/${TMP_BKG}"
        FULLREQUEST="${REQ_BASE}?${REQ_OUT}&${REQ_FILE}&${REQ_BKG}"
    else
        FULLREQUEST="${REQ_BASE}?${REQ_OUT}&${REQ_FILE}"
    fi

    echo "[STATUS] Requesting: ${FULLREQUEST}"
    if [[ ! -s ${RESULTDIR}/${file_new} ]]; then
        wget -O ${RESULTDIR}/${file_new} --waitretry=0 --tries=4 --timeout=180 "${FULLREQUEST}"
    fi
    # curl --output ${RESULTDIR}/${file_new} --retry-delay 10 --retry 4 --connect-timeout 180 "${FULLREQUEST}"

    # Compress output file.
    if [[ -s ${RESULTDIR}/${file_new} ]]; then
        gzip -f ${RESULTDIR}/${file_new}
        zcat ${RESULTDIR}/${file_new}.gz | awk -F"\t" '$1 !~ /^#/{printf $1; for(i=2; i<=22; i++){printf "\t%s",$i }; printf "\n"}' | gzip -c > ${RESULTDIR}/${file_new}_noregions.gz
        ls -sh ${RESULTDIR}/${file_new}.gz
    else
        echo "[STATUS] Failed - no lines in output."
    fi
else
    echo "[STATUS] Output file already exists at ${RESULTDIR}/${file_new}.gz"
fi

rm -rf ${PUBROOT}/${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished GREAT submission in $runtime seconds."
