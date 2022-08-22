#!/bin/bash
# --------------------------------------------------------
# Convert bedgraph to bw according to public data release:
# Put all of these into a directory on compbio_ce
# --------------------------------------------------------
# Run as: 
# qsub -cwd -t 1-17938 -l h_vmem=4G -l h_rt=00:10:00 -N conv_wg_bdg -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/make_bw_for_release.sh"
# qsub -cwd -t 1-3 -l h_vmem=5G -l h_rt=00:10:00 -N conv_wg_bdg -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/make_bw_for_release.sh"
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

STOREDIR=/broad/compbio_ce/cboix
CEBWDIR=${STOREDIR}/EPIMAP_ANALYSIS/bwtracks

# NOTE: Put each into its own mark directory:
if [[ "$( echo ${FILE} | grep "impute" - -c )" == "1" ]]; then
    bdgdir=${IMPUTED_DIR}
    BWDIR=${PUBLIC_BWDIR}/imputed/${MARK}
    BDGDIR=${PUBLIC_BDGDIR}/imputed/${MARK}
    OUTDIR=${CEBWDIR}/imputed/${MARK}
    FINALBW=${OUTDIR}/imputed_${SAMPLE}_${MARK}.bigWig
else
    bdgdir=${CONVERTED_DIR}
    BWDIR=${PUBLIC_BWDIR}/observed/${MARK}
    BDGDIR=${PUBLIC_BDGDIR}/observed/${MARK}
    OUTDIR=${CEBWDIR}/observed/${MARK}
    FINALBW=${OUTDIR}/observed_${SAMPLE}_${MARK}.bigWig
fi
mkdir -p ${BWDIR} ${OUTDIR}

if [[ ! -s $FINALBW ]]; then
    tmpwig=${TMP_DIR}/${FILE}_full.tmp.wig
    rm $tmpwig
    while read chr size; do 
        echo ${chr}
        fname=${chr}_${FILE}.wig.gz
        cfile=${bdgdir}/$fname
        zcat $cfile | awk 'NR > 1' >> $tmpwig
        # Copy the BW file to the public bedgraph directory as well.
        cp $cfile $BDGDIR/$name
    done < ${CHROMSIZES_noY}

    # Clip when converting because some 25 bp intervals 
    # in fixedStep may run over the tail.
    wigToBigWig -clip $tmpwig ${CHROMSIZES_noY} $FINALBW
fi

rm -rf $TMP_DIR

end=`date +%s`
runtime=$((end-start))
echo "Finished wig to bw aggregation and conversion in $runtime seconds."
