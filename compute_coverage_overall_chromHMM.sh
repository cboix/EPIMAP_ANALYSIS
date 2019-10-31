#!/bin/bash
# ----------------------------------------
# Compute the genome coverage of one state
# across all epigenomes
# ----------------------------------------
NSTATES=18
MNAME="observed_aux_18_on_mixed_impobs_QCUT"
if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [OPTIONS] 
    Compute the genome coverage of one state across all epigenomes
    -c     Chromosome [required]
    -s     State number (if -1, then all minus Quiescent) [required]
    -n     Number of states model (15, 18, 25) [default: $NSTATES]
    -m     MODEL to use (overrides -n) (OPTIONAL)"
    exit 1
fi

while getopts c:s:n:m: o
do      case "$o" in
    c)		CHROM="$OPTARG";;
    s)		STATE="$OPTARG";;
    n)		NSTATES="$OPTARG";;
    m)		MNAME="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

echo "[STATUS] OPTIONS set to:"
echo "[STATUS] state number: $STATE"
echo "[STATUS] nstates: $NSTATES model $MNAME"

# 0. Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh ${NSTATES}
if [[ "$MNAME" != "" ]];then
    MODELNAME=$MNAME
fi
FILEDIR=${CALLDIR}/$MODELNAME
TMP_DIR=${TMP}/extract_region_${RANDOM}
mkdir -p ${TMP_DIR}
cd $FILEDIR


# 1. Intersect or union:
rm -rf ${TMP_DIR}/tmp_aggregate.bed
while read file; do 
    echo $file
    if [[ "$STATE" == "-1" ]]; then
        # zcat $file | grep "E$NSTATES" | awk -vOFS="\t" '{print $1,$2,$3}' > ${TMP_DIR}/tmp_current.bed
        # 1-4 and 7-11# EXCLUDE 5,6,12-18
        zcat $file |  awk -vOFS="\t" '$4 ~ /E5|E6|E12|E13|E14|E15|E16|E17|E18/{print $1,$2,$3}' | sort -k1,1V -k2,2n | bedtools merge -i - > ${TMP_DIR}/tmp_current.bed
        zcat $file |  awk -vOFS="\t" '$4 ~ /E5|E6|E12|E13|E14|E15|E16|E17|E18/{print $1,$2,$3,$4}' | sort -k1,1V -k2,2n | bedtools merge -i - > ${TMP_DIR}/tmp_current.bed
        if [[ ! -s ${TMP_DIR}/tmp_aggregate.bed ]]; then
            cat ${TMP_DIR}/tmp_current.bed > ${TMP_DIR}/tmp_aggregate.bed
        else
            bedtools intersect -a ${TMP_DIR}/tmp_current.bed -b ${TMP_DIR}/tmp_aggregate.bed > ${TMP_DIR}/tmp_intersect.bed
            sort -k1,1V -k2,2n ${TMP_DIR}/tmp_intersect.bed | bedtools merge -i - > ${TMP_DIR}/tmp_aggregate.bed
        fi
    else 
    fi
done < <( ls ${FILEDIR}/BSS*segments.bed.gz )

# 2. Merge/condense:
sort -k1,1V -k2,2n ${TMP_DIR}/tmp_aggregate.bed | bedtools merge -i - > ${FILEDIR}/aggregated_E${STATE}.bed

# 3. Count number of bp:
awk '{a=a+$3-$2}END{print a}' ${FILEDIR}/aggregated_E${STATE}.bed 

# Clean up directories:
if [[ -d ${TMP_DIR} ]]; then 
    rm -rf ${TMP_DIR}
fi

end=`date +%s`
runtime=$((end-start))
echo "Finished in $runtime seconds."
