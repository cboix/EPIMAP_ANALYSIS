#!/bin/bash
# ---------------------------------------
# Call narrow and gapped peaks with MACS2
# Using imputed log10pval signal tracks
# ---------------------------------------
if [[ $# -lt 1 ]]; then
    echo "USAGE: $(basename $0) [INFOFILE] (optional [OUTDIR] [INDIR] [TASK])" >&2
    echo '  [INFOFILE]: sample/mark table to use' >&2
    echo '  [OUTDIR]: (OPTIONAL) directory to put peaks' >&2
    echo '  [INDIR]: (OPTIONAL) directory where bw are located' >&2
    echo '  [TASK]: (OPTIONAL) sample/mark combination to run (line number from info table)' >&2
    exit 1
fi
start=`date +%s`
hostname -f
INFOFILE=$1

if [[ $# -gt 1 ]]; then
    OUTDIR=$2
else
    OUTDIR=${PKIMP_DIR}
fi

if [[ $# -gt 2 ]]; then
    INDIR=$3
else
    INDIR=${IMPUTED_DIR}
fi

if [[ $# -gt 3 ]]; then
    TASK=$4
else
    TASK=${SGE_TASK_ID}
fi

TMP_DIR=${TMP}/callImpPeaks_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR} ${OUTDIR}

# Get sample/mark from table:
SAMPLE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
MARK=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $2}' )
prefix=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $3}' )
if [[ "${prefix}" == "" ]]; then
    prefix=impute_${SAMPLE}_${MARK}
fi
echo "[STATUS] Impute peaks on $SAMPLE and ${MARK}. File prefix is ${prefix}."

# ---------------------------------------------------------
# 1. Merge bedgraphs and turn them into bdg track (TMP_DIR)
# ---------------------------------------------------------
TMPBDG=$TMP_DIR/${prefix}.bdg
echo "track type=bedGraph name=${prefix}" > ${TMPBDG}
echo "[STATUS] Creating temporary bedgraph from imputed bigwig format."
while read chr size; do
    echo $chr
    chrfile=${INDIR}/${chr}_${prefix}.wig.gz
    zcat ${chrfile} | awk -vOFS="\t" -v chr=${chr} 'BEGIN{a=1}NR>2{print chr,a, a+24,$1; a=a+25}' >> ${TMPBDG}
done < <( sort -k1V $CHROMSIZES_noY )

# ------------------------------------
# 2. Use bedgraph track to call peaks:
# 
# TODO: Set read peak size differently by mark: 
# based on average of qc values of observed data.
# ------------------------------------
# MACS2 options:
set +u; source activate macs2_env; PYTHONPATH=/broad/compbio/cboix/software/miniconda2/lib/python2.7/site-packages; set -u
MAXGAP=100

# Get median fragment length for this type of mark:
MINLEN=$( grep $MARK $DBDIR/qc_reduced.tsv | cut -f4 | sort -n | awk '$1 >=40 && $1 <= 500' | awk '{a[i++]=$1;}END{x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}' ) 

# a. Narrow peaks:
OUTFILE=${OUTDIR}/${prefix}_peaks.narrowPeak
if [[ ! -s ${OUTFILE}.gz ]] || 
    [[ "$( zcat ${OUTFILE}.gz | wc -l | awk '{print $1}')" -lt 2 ]]; then
    echo "[STATUS] Calling narrow peaks from bedgraph for ${prefix}"
    macs2 bdgpeakcall -i ${TMPBDG} -c 2 -l $MINLEN -g $MAXGAP -o ${prefix}_peaks.narrowPeak --outdir ${OUTDIR}
    gzip -f ${OUTFILE}
    echo "[STATUS] Called $( zcat $OUTFILE | wc -l ) narrow peaks for ${prefix}."
fi

# b. Broad peaks:
OUTFILE=${OUTDIR}/${prefix}_peaks.broadPeak
if [[ ! -s ${OUTFILE}.gz ]] ||
    [[ "$( zcat ${OUTFILE}.gz | wc -l | awk '{print $1}')" -lt 2 ]]; then
    echo "[STATUS] Calling broad peaks from bedgraph for ${prefix}"
    macs2 bdgbroadcall -i ${TMPBDG} -c 2 -l $MINLEN -g $MAXGAP -o ${prefix}_peaks.broadPeak --outdir ${OUTDIR}
    gzip -f ${OUTFILE}
    echo "[STATUS] Called $( zcat $OUTFILE | wc -l ) broad peaks for ${prefix}."
fi

# c. Narrow peaks with much smaller gap (10) - ostensibly no merge:
OUTFILE=${OUTDIR}/${prefix}_peaks.nomerge.narrowPeak
if [[ ! -s ${OUTFILE}.gz ]] ||
    [[ "$( zcat ${OUTFILE}.gz | wc -l | awk '{print $1}')" -lt 2 ]]; then
    echo "[STATUS] Calling narrow peaks with no merge from bedgraph for ${prefix}"
    macs2 bdgpeakcall -i ${TMPBDG} -c 2 -l $MINLEN -g 10 -o ${prefix}_peaks.nomerge.narrowPeak --outdir ${OUTDIR}
    gzip -f ${OUTFILE}
    echo "[STATUS] Called $( zcat $OUTFILE | wc -l ) non-merged narrow peaks for ${prefix}."
fi
# source deactivate > /dev/null

# -----------------------------------------
# 3. Refine peaks in order to create index.
# -----------------------------------------
# Intersect bedgraph with non-merged peaks:
SUMMFILE=${OUTDIR}/${prefix}_peaks.summit.narrowPeak.gz
PKFILE=${OUTDIR}/${prefix}_peaks.nomerge.narrowPeak.gz 
if [[ ! -s $SUMMFILE ]] || [[ $PKFILE -nt $SUMMFILE ]] || 
    [[ "$( zcat ${SUMMFILE} | wc -l | awk '{print $1}')" -lt 2 ]]; then
    zcat ${PKFILE} | awk -F"\t" -vOFS="\t" 'NR > 1{print $1,$2,$3,$4,$5}' | sort -k1,1V -k2,2n > ${TMP_DIR}/query.bed
    echo "[STATUS] Sorting bdg"
    sort-bed $TMPBDG > ${TMP_DIR}/sorted.bdg

    echo "[STATUS] Extracting relevant bdg areas"
    sort-bed ${TMP_DIR}/query.bed | bedextract ${TMP_DIR}/sorted.bdg - | bedtools intersect -a - -b ${TMP_DIR}/query.bed -wb | cut -f1,2,3,4,8 > ${TMP_DIR}/intersect.bed

    # Take max + its midpt:
    awk -F"\t" -vOFS="\t" 'BEGIN{pk=""; max=0; chr=""; mid=""}{if($5 == pk){if($4>max){max=$4; mid=($2 + $3)/2;}} else {if(pk!=""){print chr,mid,max,pk}; chr=$1; max=$4; mid=($2 + $3)/2; pk=$5; } }END{print chr,mid,max,pk}' ${TMP_DIR}/intersect.bed | sort -k4 > ${TMP_DIR}/proc_intersect.bed

    zcat $PKFILE | sort -k4 | join -1 4 -2 4 - ${TMP_DIR}/proc_intersect.bed | awk -vOFS="\t" '{print $2,$3,$4,$1,$13,$12}' | sort -k4V | gzip -c > $SUMMFILE

    echo "[STATUS] Called summits for $( zcat $SUMMFILE | wc -l ) peaks for $prefix"
else
    echo "[STATUS] Summits already exist for $prefix"
fi

rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished calling peaks from bw sucessfully in $runtime seconds."
