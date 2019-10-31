#!/bin/bash
# SGE Array to process ChIP-seq epitopes
# Note should use this to process large batch epitopes in parallel
export CELL=$( sed "${SGE_TASK_ID}q;d" ${EPDIR}/cell_types )
export SCELL=$CELL # saved as it may get changed.
export NREADS=30000000 # Subsampling depth

# ====================================
# STEP1 - Download + Process each one:
# ====================================
RERUNSTEP2=0 # If we update step 1
IFS=$'\t'
while read link rep assembly processing cell file
do
    # Get attributes:
    id=$(echo "$link" | awk '{sub(".*download/","",$0); sub("\\..*$","",$0); print $0}' )
    if [[ "$processing" == "unfiltered alignments" ]] ; then
        processing="unfiltered_alignments"
    fi
    echo "$id $link $assembly $cell $processing"

    # Download and pre-process to tagAlign:
    STEP1_FILE="${TADIR}/FINAL_${id}_${cell}.tagAlign.gz"
    if LC_ALL=C gzip -l ${STEP1_FILE} | awk 'NR==2 {exit($2!=0)}'
    then
        echo "- STEP1 Preprocess $id $cell $EPITOPE"
        echo "- Making ${STEP1_FILE}"
        source $BINDIR/BAMtoTA_download_preprocess.sh $id $cell $link $assembly $processing
        # if updates happen, run step 2 again:
        RERUNSTEP2=1
    fi
done < <( awk -vFS="\t" -v cell=$CELL '$5 == cell' $FILEINFO ) 

# ==========================================
# STEP2 - Pool same cell type + perform SCCA
# Output a bedgraph for subsampled + not
# ==========================================
# Check if pooled contains all available reads: 
STEP2_FILE="${TADIR}/FINAL_${EPITOPE}_${CELL}.tagAlign.gz"
NTOT=$( cat ${QCDIR}/FINAL_${EPITOPE}_${CELL}.numreads) 
NSUM=$( cat ${QCDIR}/ENCFF*${CELL}_*.numreads | awk 'BEGIN{a=0}{a=a+$1}END{print a}' ) 
[[ "$NTOT" != "$NSUM" ]] && RERUNSTEP2=1
if [[ "$( LC_ALL=C gzip -l ${STEP2_FILE} | awk 'NR==2{print $2}')" == "0" || "$RERUNSTEP2" == "1" ]]
then
    echo "- STEP2 Pool all files for ${CELL} in ${EPITOPE}"
    echo "- Making $STEP2_FILE"
    source $BINDIR/Pool_subsample_SCCA.sh $CELL $NREADS
    echo "STEP2: SUBSAMPLED to $( zcat ${STEP2_FILE} | wc -l ) reads"
fi

# Stop if WCE - doesn't need to go further
if [[ "${EPITOPE}" == "WCE" ]] 
then
    echo "DONE: WCE only processed and pooled, no peak calling."
    exit 1
fi

# ====================================
# STEP3 - Peak + broadPeak calling:
# ====================================
RERUNSTEP3=$RERUNSTEP2
controlstub="FINAL_WCE_${CELL}" 
chipstub="FINAL_${EPITOPE}_${CELL}.sub"
CPREF=${chipstub}_VS_${controlstub}
outPref="${PKDIR}/${CPREF}" 
STEP3_FILE="${outPref}_peaks.narrowPeak.gz"
# Check if WCE FILE UPDATED:
CONTROL_FILE=${CHPDIR}/files/WCE/tagAlign/${controlstub}.tagAlign.gz

[[ "${CONTROL_FILE}" -nt "${STEP3_FILE}" ]] && RERUNSTEP3=1
if [[ "$( LC_ALL=C gzip -l ${STEP3_FILE} | awk 'NR==2 {print $2}' )" == "0" || "${RERUNSTEP3}" == "1" ]] 
then
    echo "- STEP3 Call narrow and broad peaks for ${CELL} in ${EPITOPE}"
    echo "- Making $STEP3_FILE"
    source $BINDIR/Call_peaks.sh $CELL
fi

# ====================================
# STEP4 - Signal Track generation:
# ====================================
outPref="${BDGDIR}/${CPREF}" 
STEP4_FILE=${outPref}.pval.signal.bedgraph.gz
if LC_ALL=C gzip -l ${STEP4_FILE} | awk 'NR==2 {exit($2!=0)}'
then
    echo "- STEP4 Generate signal tracks for ${CELL} in ${EPITOPE}"
    echo "- Making $STEP4_FILE"
    source $BINDIR/Generate_signal_tracks.sh $CELL
fi

# =============================================
# STEP4B - Convert pval files for ChromImpute: 
# =============================================
STEP4B_FILE=${CONVERTED_DATADIR}/${EPITOPE}/chr20_${CPREF}.pval.signal.wig.gz
if LC_ALL=C gzip -l ${STEP4B_FILE} | awk 'NR==2 {exit($2!=0)}'
then
    echo "Convert signal into 25bp average tracks for ChromImpute"
    echo "${cell}\t${EPITOPE}\t${chipstub}_VS_${controlstub}.pval.signal.bedgraph.gz" > ${outPref}_tab.txt
    java -mx4000M -jar ${CHROMIMPUTE} Convert ${BDGDIR} ${outPref}_tab.txt ${CHROMSIZES} ${CONVERTED_DATADIR}/${EPITOPE}
    rm ${outPref}_tab.txt
fi

# ===========================
# STEP5 - Pseudo-reps and IDR
# ===========================
