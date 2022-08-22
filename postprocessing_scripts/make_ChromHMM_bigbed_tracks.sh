#!/bin/bash
# ---------------------------------------------------
# Make bigbed files for the hg19/hg38 ChromHMM files:
# ---------------------------------------------------
# Directories/Config:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh
start=`date +%s`
hostname -f

# Arguments:
TMP_DIR=${TMP}/chromhmm_to_bigbed_${RANDOM}
mkdir -p ${TMP_DIR} $STOREDIR

COLTSV=${DBDIR}/CHMM_${MODELNAME}_colors.tsv
TMPCOL=${TMP_DIR}/CHMM_colors.tsv
awk -vOFS="\t" -F"\t" '{print "E"$1, $1"_"$2,$3}' $COLTSV | sort -k1 > $TMPCOL

for SGE_TASK_ID in `seq 1 833`; do
    echo ${SGE_TASK_ID}
    datasets=(observed_aux_18_hg19 observed_aux_18_hg38)
    for dset in ${datasets[@]}; do
        echo $dset
        CALLDIR=${PUBLIC_CMDIR}/${dset}/CALLS
        BBDIR=${CALLDIR}/bigBed/
        mkdir -p ${BBDIR} 

        cd $CALLDIR
        sfile=$( ls BSS*segments.bed.gz | sed "${SGE_TASK_ID}q;d" )
        sprefix=${sfile%%.bed.gz}
        echo "[STATUS] Converting ${sfile}"
        BBFILE=${BBDIR}/${sprefix}.bb

        # Set chrom sizes file:
        if [[ "$dset" == "observed_aux_18_hg19" ]]; then
            GENSIZE=$ANNDIR/hg19.chrom.sizes_noY
        else 
            GENSIZE=$ANNDIR/hg38.chrom.sizes_noY
        fi

        if [[ ! -s $BBFILE ]]; then
            # Convert to BED9
            # zcat $sfile | awk '$1 ~ /^chr[0-9XY]*$/' | sort -k4 | join -1 4 -2 1 - $TMPCOL > ${TMP_DIR}/tmp.bed
            # awk -vOFS="\t" '{print $2,$3,$4, $5,"0",".",$3,$4,$6}' ${TMP_DIR}/tmp.bed | sort -k1,1 -k2,2n > ${TMP_DIR}/tmp.srt.bed
            # Already BED9, skip above:
            zcat $sfile | sort -k1,1 -k2,2n > ${TMP_DIR}/tmp.srt.bed
            bedToBigBed ${TMP_DIR}/tmp.srt.bed ${GENSIZE} ${BBFILE}
        fi
    done
done

rm -rf $TMP_DIR

end=`date +%s`
runtime=$((end-start))
echo "Finished wig to bw aggregation and conversion in $runtime seconds."
