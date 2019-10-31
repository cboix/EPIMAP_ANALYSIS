#!/bin/bash
# SGE Array to process ChIP-seq epitopes
# Note should process small batch epitopes - large ones should be parallelized.
export EPITOPE=$1
export FILEINFO=${CHPLNK}/${EPITOPE}_bam.csv
export EPDIR=${CHPDIR}/files/${EPITOPE}
export TADIR=${EPDIR}/tagAlign
export PKDIR=${EPDIR}/peaks
export BDGDIR=${EPDIR}/bedgraph
export QCDIR=${EPDIR}/qc
mkdir -p $EPDIR $TADIR $QCDIR $BDGDIR $PKDIR
echo "Processing epitope: $EPITOPE - has $( wc -l $FILEINFO | awk '{print $1 - 1}' ) files"

# ====================================
# STEP1 - Download + Process each one:
# ====================================
IFS=$'\t'
while read link rep assembly processing cell file
do
    # Get attributes:
    id=$(echo $link | awk '{sub(".*download/","",$0); sub("\\..*$","",$0); print $0}' )
    if [[ "$processing" == "unfiltered alignments" ]] ; then
        processing="unfiltered_alignments"
    fi
    echo "$id $link $assembly $cell $processing"

    STEP1_FILE="${TADIR}/FINAL_${id}_${cell}.tagAlign.gz"
    if LC_ALL=C gzip -l ${STEP1_FILE} | awk 'NR==2 {exit($2!=0)}'
    then
        echo "- STEP1 Preprocess $id $cell $EPITOPE"
        echo "- Making ${STEP1_FILE}"
        source $BINDIR/BAMtoTA_download_preprocess.sh $id $cell $link $assembly $processing
    fi
done < <(awk 'NR>1' $FILEINFO) # skip header

# ==========================================
# STEP2 - Pool same cell type + perform SCCA
# Output a bedgraph for subsampled + not
# ==========================================
echo "STEP 2 for ${epitope}"
awk -vFS="\t" 'NR>1{print $5}' $FILEINFO | sort -u > ${EPDIR}/cell_types
IFS=$'\t'
while read cell
do 
    STEP2_FILE="${TADIR}/FINAL_${EPITOPE}_${cell}.tagAlign.gz"
    if LC_ALL=C gzip -l ${STEP2_FILE} | awk 'NR==2 {exit($2!=0)}'
    then
        echo "- STEP2 Pool all files for ${cell} in ${EPITOPE}"
        echo "- Making $STEP2_FILE"
        source $BINDIR/Pool_subsample_SCCA.sh $cell
    fi
done < ${EPDIR}/cell_types 

# ====================================
# STEP3 - Peak calling:
# Perform IDR and generate peak calls
# Generate binned data in CHMM format.
# ====================================

