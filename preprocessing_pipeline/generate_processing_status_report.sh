#!/bin/bash 
# NOTE: Run vars from main_EPIMAP_ANALYSIS.sh first.
# ====================================
# Create a report of processing status
# for debugging and tracking processing.
# ====================================
# TODO GOAL: Create Images + collate into pdf report

# ===========================
# PROCESSING INFO TO INCLUDE: 
# ===========================
# Tables of CELL TYPES vs MARKS
# For all and for just histone marks.
# 1. Number of replicates
# 2. Number of reads total
# 3. Number of reads in each replicate
# 4. Processing status

# ==========================
# 1. Generate report matrix: 
# ==========================
TODAY=$(date +%Y%m%d-%H%M)
REPMAT=${ANNDIR}/processing_status_matrix_${TODAY}.tsv

echo "Epitope CellType TotalReplicates CellRepl RepsProcessed PoolSCCA PeakCalls SignalTracks Depth SubsampleDepth Fraglen FraglenSubsample" > $REPMAT

# ChIP-SEQ:
while read EPITOPE
do
    echo "Adding $EPITOPE to the report:"

    # Directories:
    if [[ "$EPITOPE" == "DNase-seq" ]]
    then
        FILEINFO=${DHSDIR}/file_links/all_submitted_released/DNase-seq_bam.tsv
        EPDIR=${DHSDIR}/files/DNase-seq
        rm -f ${QCDIR}/DNase-seq_read_totals.tsv
    else
        FILEINFO=${CHPLNK}/${EPITOPE}_bam.csv
        EPDIR=${CHPDIR}/files/${EPITOPE}
    fi
    TADIR=${EPDIR}/tagAlign
    PKDIR=${EPDIR}/peaks
    BDGDIR=${EPDIR}/bedgraph
    QCDIR=${EPDIR}/qc
    mkdir -p $EPDIR $TADIR $QCDIR $BDGDIR $PKDIR

    # Total number of files:
    NTOT=$( wc -l $FILEINFO | awk '{print $1 - 1}' )
    # Cell types and replicates: 
    awk -vFS="\t" 'NR>1{print $5}' $FILEINFO | sort -u > ${EPDIR}/cell_types
    while read cell 
    do
        # Number of replicates: 
        NREP=$( awk -vFS="\t" -v cell=$cell '$5 == cell' $FILEINFO | wc -l )
        CPREF=FINAL_${EPITOPE}_${cell}
        if [[ "$EPITOPE" == "DNase-seq" ]] 
        then
            OPREF="${CPREF}.sub_VS_Uniform_BKG_CONTROL_36_50000000" 
        else
            OPREF="${CPREF}.sub_VS_FINAL_WCE_${cell}" 
        fi

        # Determine processing step:
        STEP2_FILE="${TADIR}/${CPREF}.sub.tagAlign.gz"
        STEP3_FILE="${PKDIR}/${OPREF}.narrowPeak.gz"
        STEP4_FILE="${CONVERTED_DATADIR}/chr1_${OPREF}.pval.signal.bedgraph.gz.wig.gz"

        # How many replicates processed:
        STEP1_DONE=$( ls -s ${TADIR}/FINAL_ENCFF*_${cell}.tagAlign.gz 2>/dev/null | awk 'BEGIN{a=0}{b=($1==34?0:1); a=a+b}END{print a}' ) 
        STEP2_DONE=$( ls -s ${STEP2_FILE} 2>/dev/null | awk 'BEGIN{a=0}{a=($1==34?0:1)}END{print a}' )
        STEP3_DONE=$( ls -s ${STEP3_FILE} 2>/dev/null | awk 'BEGIN{a=0}{a=($1==34?0:1)}END{print a}' )
        STEP4_DONE=$( ls -s ${STEP4_FILE} 2>/dev/null | awk 'BEGIN{a=0}{a=($1==34?0:1)}END{print a}' )

        # Total and subsampled depth and est. fragment length (from stats)
        if [[ "${STEP2_DONE}" == "1" ]] 
        then
            STEP1_DONE=1  # May have removed files.
            NDEPTH=$( awk '{printf "%f", $1 / 1000000}' ${QCDIR}/${CPREF}.numreads )
            NSUB=$( awk '{printf "%f", $1 / 1000000}' ${QCDIR}/${CPREF}.sub.numreads )
            # FRAGTOT=$( awk '{print $3}' ${QCDIR}/${CPREF}.cc.qc )
            FRAGLEN=$( awk '{print $3}' ${QCDIR}/${CPREF}.sub.cc.qc )
            FRAGTOT="$FRAGLEN"
        else 
            NDEPTH=0; NSUB=0; FRAGTOT=0; FRAGLEN=0;
        fi

        # Output all: 
        echo "$EPITOPE $cell $NTOT $NREP ${STEP1_DONE} ${STEP2_DONE} ${STEP3_DONE} ${STEP4_DONE} ${NDEPTH} ${NSUB} ${FRAGTOT} ${FRAGLEN}" >> $REPMAT
    done < ${EPDIR}/cell_types
done < <( ls $CHPLNK | awk '$1 ~ /csv/{gsub("_bam.csv","",$0); print $0}' | awk '{print $0}END{print "DNase-seq"}' )

# =========================================================
# 2. Make report tables and figures from processing status:
# =========================================================
# R processing: 
R --slave -f $BINDIR/plot_processing_status.R

# =============================
# 3. VALIDATION INFO TO INCLUDE
# =============================
# Distance metrics different datasets
# Trees between different samples

