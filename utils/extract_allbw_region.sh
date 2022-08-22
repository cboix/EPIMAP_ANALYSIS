#!/bin/bash
# ----------------------------------------
# Extract a region across all BW:
# observed + imputed, in all marks.
# Turn region into set of legible matrices
# ----------------------------------------
if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [OPTIONS] 
    -c     Chromosome [required]
    -a     Start coordinate [required]
    -b     End coordinate [required]
    -o     Output prefix [required]"
    exit 1
fi

CHROM="chr15"
# CSTART="26927000"
# CEND="26933000"
CSTART="27827000"
CEND="27835000"

OUTPREF=${TMP}/testchunk_${CHROM}_${CSTART}_${CEND}
while getopts c:a:b:o: o
do      case "$o" in
    c)		CHROM="$OPTARG";;
    a)		CSTART="$OPTARG";;
    b)		CEND="$OPTARG";;
    o)		OUTPREF="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

echo "[STATUS] OPTIONS set to:"
echo "[STATUS] Region: ${CHROM}:${CSTART}-$CEND"
echo "[STATUS] Prefix: ${OUTPREF}"

# 0. Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 

TMP_DIR=${TMP}/extract_region_bw_${RANDOM}
mkdir -p ${TMP_DIR}

# -----------------------------------
# 1. Collate all corresponding files:
# TODO: PARALLELIZE BY EPITOPE?
# -----------------------------------
# Get bin start and end (1-indexed) + add two for header:
BSTART=$( echo "$CSTART" | awk '{a=$1 % 25; print ($1 - a) / 25 + 1 + 2}' )
BEND=$( echo "$CEND" | awk '{a=$1 % 25; print ($1 - a) / 25 + 2}' )
BLEN=$( echo "$BSTART $BEND" | awk '{print $2 - $1 + 1}' )
echo "[STATUS] Pulling row #s: $BSTART to $BEND, with length $BLEN"

AGTAB=${DBDIR}/best_worst_impobs_agreement_011119.tsv
rm -f ${TMP_DIR}/bw_files_*.tsv ${TMP_DIR}/concat_tmp_*.mtx
while read mark; do
    echo "----------- $mark -----------"
    while read cell epitope fpref; do
        FILELIST=${TMP_DIR}/bw_files_${epitope}.tsv
        CONCATFILE=${TMP_DIR}/concat_tmp_${epitope}.mtx
        printf "."
        if [[ "$fpref" == "impute_${cell}_${epitope}" ]];then 
            BWFILE=${IMPUTED_DIR}/${CHROM}_${fpref}.wig.gz
            data="i"
        else
            # Check quality of data:
            stats=$( grep -P "${cell}\t${epitope}\t" $AGTAB | awk '{print $3}' )
            if [[ "$stats" == "LOW" ]]; then 
                echo "${fpref} is low."
                BWFILE=${CONVERTED_DIR}/EMPTYFILE.wig.gz
                data="o"
            else
                BWFILE=${CONVERTED_DIR}/${CHROM}_${fpref}.wig.gz
                data="o"
            fi
        fi
        if [[ -s ${BWFILE} ]]; then
            echo "$cell $epitope $BWFILE" >> $FILELIST
            zcat $BWFILE | sed -n "${BSTART},${BEND}p;${BEND}q" | awk -vOFS="\t" -v pref="${cell}\t${data}" '{print NR,$1,pref }' >> ${CONCATFILE}
        fi
    done < <( grep $mark ${ALL_TRACKS_TAB} )
done < <( echo "DNase-seq\nH3K27ac\nH3K4me1\nH3K4me3\nH3K36me3\nH3K9me3" )


# Turn files into 
MARK=DNase-seq
FILELIST=${TMP_DIR}/bw_files_${MARK}.tsv
CONCATFILE=${TMP_DIR}/concat_tmp_${MARK}.mtx




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

# 3. Plot region with R:
# Arguments: dataframe and a tagline (date)
echo "[STATUS] Plotting imputed regions"
source activate mv_env > /dev/null
R --slave -f ${BINDIR}/plot_chromHMM_region.R --args ${TMP_DIR}/imputed_region_calls.bed imputed_${MODELNAME}_${TODAY}
source deactivate > /dev/null

# TODO repeat with ROADMAP DATA:

# Clean up directories:
if [[ -d ${TMP_DIR} ]]; then 
    rm -rf ${TMP_DIR}
fi

end=`date +%s`
runtime=$((end-start))
echo "Finished enrichment pipeline in $runtime seconds."
