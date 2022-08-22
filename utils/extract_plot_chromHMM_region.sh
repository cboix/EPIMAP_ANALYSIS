#!/bin/bash
# -----------------------------------
# Extract and plot a chromHMM region:
# 1) Imputed
# 2) Roadmap
# 3) Comparison of same cells in both
# -----------------------------------
# NSTATES=15 
# MNAME=""
NSTATES=18
MNAME="observed_aux_18_on_mixed_impobs_QCUT"
RDMAP=1
if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [OPTIONS] 
    -c     Chromosome [required]
    -a     Start coordinate [required]
    -b     End coordinate [required]
    -r     Add original roadmap to plot (OPTIONAL 0/1) [default: $RDMAP]
    -n     Number of states model (15, 18, 25) [default: $NSTATES]
    -m     MODEL to use (overrides -n) (OPTIONAL)"
    exit 1
fi

while getopts c:a:b:n:m:r: o
do      case "$o" in
    c)		CHROM="$OPTARG";;
    a)		CSTART="$OPTARG";;
    b)		CEND="$OPTARG";;
    n)		NSTATES="$OPTARG";;
    m)		MNAME="$OPTARG";;
    r)      RDMAP="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

echo "[STATUS] OPTIONS set to:"
echo "[STATUS] chr $CHROM start $CSTART end $CEND"
echo "[STATUS] nstates $NSTATES model $MNAME"

# 0. Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh ${NSTATES}
if [[ "$MNAME" != "" ]];then
    MODELNAME=$MNAME
fi
FILEDIR=${CALLDIR}/$MODELNAME
TMP_DIR=${TMP}/extract_region_${RANDOM}
mkdir -p ${TMP_DIR}
cd $FILEDIR

# 1. Extract calls using bedtools intersect: 
BEDFILE=${TMP_DIR}/original_region.bed
echo -e "${CHROM}\t${CSTART}\t${CEND}" > $BEDFILE
cat $BEDFILE

sort -k1V ${CHROMSIZES} >  ${TMP_DIR}/genome.tmp
intersect_cmd="bedtools intersect -a ${BEDFILE} -b $( ls *segments.bed.gz | awk '{printf $0" "}' ) -sorted -wb -filenames -g ${TMP_DIR}/genome.tmp > ${TMP_DIR}/intersect_region.bed"

echo "[STATUS] Running bedtools intersect to identify regions"
echo "$intersect_cmd"
bash -c "$intersect_cmd"

# 2. Clean up intersect file:
awk -vOFS="\t" '{sub("_[1-3][0-9]_CALLS_segments.bed.gz","",$4); print $1,$2,$3,$4,$8}' ${TMP_DIR}/intersect_region.bed > ${TMP_DIR}/imputed_region_calls.bed

if [[ "$RDMAP" == "1" ]]; then
    cd $RDCHMMDIR/n$NSTATES

    intersect_cmd="bedtools intersect -a ${BEDFILE} -b $( ls *segments.bed.gz | awk '{printf $0" "}' ) -sorted -wb -filenames -g ${TMP_DIR}/genome.tmp > ${TMP_DIR}/roadmap_intersect_region.bed"

    echo "[STATUS] Running bedtools intersect to identify regions"
    echo "$intersect_cmd"
    bash -c "$intersect_cmd"

    # 2b. Clean up intersect file:
    awk -vOFS="\t" '{sub("_.*_segments.bed.gz","",$4); print $1,$2,$3,$4,$8}' ${TMP_DIR}/roadmap_intersect_region.bed > ${TMP_DIR}/roadmap_region_calls.bed
fi

# Get the genes:
# zcat $ANNDIR/gencode.v30lift37.basic.annotation.gtf.gz | awk '$0 ~ /protein_coding/ && ($3 ~ /gene|exon/)' | awk -vOFS="\t" -v chr=$CHROM -v cend=$CEND -v cst=$CSTART '$1 == chr && ($4 < cend || $5 > cst)' > tmp


# 3. Plot region with R:
# Arguments: dataframe and a tagline (date)
echo "[STATUS] Plotting imputed regions"
source activate mv_env > /dev/null
if [[ "$RDMAP" == "1" ]]; then
    R --slave -f ${BINDIR}/plot_chromHMM_region.R --args ${TMP_DIR}/imputed_region_calls.bed imputed_${MODELNAME} $TMP_DIR/roadmap_region_calls.bed
else
    R --slave -f ${BINDIR}/plot_chromHMM_region.R --args ${TMP_DIR}/imputed_region_calls.bed imputed_${MODELNAME}
fi
source deactivate > /dev/null

# TODO repeat with ROADMAP DATA:

# Clean up directories:
if [[ -d ${TMP_DIR} ]]; then 
    rm -rf ${TMP_DIR}
fi

end=`date +%s`
runtime=$((end-start))
echo "Finished enrichment pipeline in $runtime seconds."
