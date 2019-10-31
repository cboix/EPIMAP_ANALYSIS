#!/bin/bash
# --------------------------------------------------------
# Convert bedgraph to bw according to public data release:
# AND transfer to the Wang Lab data exchange server.
# --------------------------------------------------------
# Run as: 
# qsub -cwd -P compbio_lab -t 1-17938 -tc 500 -l h_vmem=4G -l h_rt=00:20:00 -N conv_wg_bdg -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/make_and_transfer_bw_for_release.sh"
# Directories/Config:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh
start=`date +%s`
hostname -f

# Arguments:
TASK=${SGE_TASK_ID}
INFOFILE=$PUBLIC_METADIR/all_released_tracks.tsv
TMP_DIR=${TMP}/bdg_bw_release_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR}

# Get sample/mark from table (available or to_impute):
SAMPLE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
MARK=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $2}' )
FILE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $3}' )
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
STATUSFILE=$BDGDIR/${SAMPLE}_${MARK}.status.txt
mkdir -p $BWDIR $BDGDIR

if [[ ! -s $STATUSFILE ]] || [[ "$( cat $STATUSFILE)" != "0" ]]; then
    tmpwig=${TMP_DIR}/${FILE}.wig
    rm $tmpwig
    while read chr size; do 
        echo ${chr}
        fname=${chr}_${FILE}.wig.gz
        cfile=${bdgdir}/$fname
        zcat $cfile | awk 'NR > 1' >> $tmpwig
        # Copy the BW file to the public bedgraph directory as well.
    done < ${CHROMSIZES_noY}

    gzip ${tmpwig}
    ls -sh ${tmpwig}.gz

    # SFTP Transfer command:
    export BATCHFILE=${TMP_DIR}/batchfile.sh
    echo "cd upload/imputation/${bdgclass}/" > $BATCHFILE
    echo "put ${tmpwig}.gz" >> $BATCHFILE
    echo "bye" >> $BATCHFILE
    cat $BATCHFILE

    # Send file over:
    export SSHPASS="share#123$"
    sshpass -e sftp -oBatchMode=no -b $BATCHFILE -oPort=22 wangshare@target.wustl.edu

    # Log the exitcode:
    echo "$?" > ${STATUSFILE}
    cat ${STATUSFILE}
fi

rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished wig compilation and send in $runtime seconds."
