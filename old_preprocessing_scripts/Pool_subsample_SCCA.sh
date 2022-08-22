#!/bin/bash -l
# =====================================================================
# Pool all replicates and subsample to 30M READS, according to Roadmap.
# We also run phantompeakqualtools in order to obtain qc stats and SCCA
# Finally, we turn the files into bedgraphs from ChromImpute & plotting
# =====================================================================
cell=$1
NREADS=30000000 # subsampling depth
if [[ $# -gt 1 ]]; then NREADS=$2; fi

# Prefixes:
CPREF=FINAL_${EPITOPE}_${cell}
QCPREFIX=${QCDIR}/${CPREF}
OFPREFIX=${TADIR}/${CPREF}
FINAL_TA_FILE="${OFPREFIX}.tagAlign.gz"
FINAL_SUB_FILE="${OFPREFIX}.sub.tagAlign.gz"

# =================================
# Pool all tagAlign files for the same cell type and epitope
# Then subsample to 30M reads
# ================================
if LC_ALL=C gzip -l ${FINAL_TA_FILE} | awk 'NR==2 {exit($2!=0)}'
then
    echo "Pool files"
    zcat ${TADIR}/FINAL_ENCFF*_${cell}.tagAlign.gz | gzip -c > ${FINAL_TA_FILE}

    # Count total number of reads in TA file!
    zcat ${FINAL_TA_FILE} | wc -l > ${QCPREFIX}.numreads
    NUMR=$(cat ${QCPREFIX}.numreads )
    echo "Have ${NUMR} reads and require ${NREADS} for subsampling."
    
    # FIXME must deal if not enough reads.
    echo "Subsample to create ${FINAL_TA_FILE}"
    zcat ${FINAL_TA_FILE} | shuf -n ${NREADS} | sort -k1,1V -k2,2g | gzip -c > ${FINAL_SUB_FILE}
fi 

# ========================
# Use SPP for SCCA and QC:
# ========================
if [[ ! -s ${QCPREFIX}.cc.qc || ! -s ${QCPREFIX}.sub.cc.qc ]] 
then
    # SCCA: 
    echo "Running SPP for SCCA on the ${CPREF} TA files."
    Rscript ${SPPDIR}/run_spp_nodups.R -c=${FINAL_TA_FILE} -filtchr=chrM -savp=${QCPREFIX}.cc.plot.pdf -out=${QCPREFIX}.cc.qc
    sed -i -r 's/,[^\t]+//g' ${QCPREFIX}.cc.qc # Keeps best length.

    # SCCA on subset:
    Rscript ${SPPDIR}/run_spp_nodups.R -c=${FINAL_SUB_FILE} -filtchr=chrM -savp=${QCPREFIX}.sub.cc.plot.pdf -out=${QCPREFIX}.sub.cc.qc
    sed -i -r 's/,[^\t]+//g' ${QCPREFIX}.sub.cc.qc

    echo "QC stats on the final tagAlign files."
    Rscript ${SPPDIR}/run_spp_nodups.R -rf -c=${FINAL_TA_FILE} -savp=${QCPREFIX}.tagAlign.pdf -out=${QCPREFIX}.stats.qc
    Rscript ${SPPDIR}/run_spp_nodups.R -rf -c=${FINAL_SUB_FILE} -savp=${QCPREFIX}.sub.tagAlign.pdf -out=${QCPREFIX}.sub.stats.qc
fi

