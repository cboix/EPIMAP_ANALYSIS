#!/bin/bash
# -------------------------------------------
# Make enhancer state list for a single model
# -------------------------------------------
# Directories/Config:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh
start=`date +%s`
hostname -f

# Arguments:
ELEMENT=ENH
MODEL=observed_aux_18_on_mixed_impobs_QCUT

# Directories:
TMP_DIR=${TMP}/make_statelist_${RANDOM}
mkdir -p ${TMP_DIR}

# Get relevant states:
CHMM_STATEMAP=${ANNDIR}/chmm_element_mapping.txt
stline=$( awk -v elem=$ELEMENT -v ns=$NUMSTATES '$1 == elem && $2 == ns' $CHMM_STATEMAP )
STATES=$( echo $stline | awk '{print $3}' )
MERGED_STATES=$( echo $stline | awk '{print $4}' )
STREG=$( echo $stline | awk '{gsub(",","|E", $3); print "E"$3}' )
STNUM=$( echo $stline | awk '{gsub(",","$|^", $3); print "^"$3"$"}' )

echo "[STATUS] Creating enh list for $ELEMENT with states: $STATES"

CDIR=${CALLDIR}/${MODEL}
ALLENH=${TMP_DIR}/all_concat.bed

# cd $CDIR
# for file in `ls BSS*segments.bed`; do
#     echo $file
#     base=$(basename $file)
#     grep -E -- "$STREG" $file | sort -k1V -k2n > ${TMP_DIR}/${base}.tmp
#     # Number regions: # TODO: Black list some afterwards?
#     awk -v base=${base%%_*} '{a=$3-$2+a}END{print a/200, base}' ${TMP_DIR}/${base}.tmp >> ${TMP_DIR}/sample_nregions.tsv
# done
# # Only chosen tracks:
# awk -vOFS="\t" '{print $2,$1}' ${TMP_DIR}/sample_nregions.tsv | sort +0 -1 > ${TMP_DIR}/samp.tsv
# sort +0 -1 $PUBLIC_METADIR/samplelist_833_final.tsv | join ${TMP_DIR}/samp.tsv - | awk -vOFS="\t" '{print $2,$1}' > ${TMP_DIR}/sample_nregions.tsv


cd $CDIR
rm ${TMP_DIR}/sample_nregions.tsv $ALLENH
while read sample; do 
    echo $sample
    ENHFILE=${CDIR}/${sample}_enhloc.bed
    if [[ ! -s $ENHFILE ]]; then
        while read chrom size; do 
            file=STATEBYLINE/${sample}_18_CALLS_PER_LINE_${chrom}_statebyline.txt.gz
            echo $file
            base=$(basename $file)
            zcat $file | grep -n -E -- "$STNUM" - | awk -vOFS="\t" -v chrom=$chrom '{split($0,a, ":"); print chrom,((a[1] - 1) * 200), (a[1] * 200),"E"a[2]}' >> ${ENHFILE}
            # sort -k1V -k2n > ${TMP_DIR}/${base}.tmp
            # Number regions: # TODO: Black list some afterwards?
        done < ${CHROMSIZES_noY}
    fi
    wc -l $ENHFILE | awk -vOFS="\t" -v samp=$sample '{print $1, samp}' >> ${TMP_DIR}/sample_nregions.tsv
    cat $ENHFILE >> ${ALLENH}
done < $PUBLIC_METADIR/samplelist_833_final.tsv 

# Stats on number of ENH 200bp bins for each sample:
wc ${TMP_DIR}/sample_nregions.tsv
sort -nr ${TMP_DIR}/sample_nregions.tsv | tail 
sort -n ${TMP_DIR}/sample_nregions.tsv | tail 

# Sort regions + make bedgraph:
SRTENH=${ALLENH%%bed}srt.bed
sort -k1V -k2n $ALLENH > $SRTENH


bedtools genomecov -i ${SRTENH} -g ${CHROMSIZES} -bg > ${TMP_DIR}/all_cov.bedgraph

cut -f4 ${TMP_DIR}/all_cov.bedgraph | sort -n | uniq -c > ${TMP_DIR}/all_hist_counts.tsv

head -n 20 ${TMP_DIR}/all_hist_counts.tsv

ls -sh ${TMP_DIR}/all_cov.bedgraph
wc ${TMP_DIR}/all_cov.bedgraph

awk '$4 > 1' ${TMP_DIR}/all_cov.bedgraph | wc
awk '$4 > 2' ${TMP_DIR}/all_cov.bedgraph | wc



# Chunk into 200 bp segments
# Choose which 200 bp segments
# Give uq id names:
# Create other auxiliary files (like for DML_DIR)


rm -rf $TMP_DIR

end=`date +%s`
runtime=$((end-start))
echo "Finished enhancer state list creation in $runtime seconds."
