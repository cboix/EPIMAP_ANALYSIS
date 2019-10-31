#!/bin/bash -l
# =====================================================================
# Call peaks in the ROADMAP/ENCODE processing pipeline (for imputation)
#   1. NarrowPeaks at pval < 0.01
#   2. BroadPeaks at pval < 0.1 containing a narrowPeak at pval < 0.01
# 
#   Adapted from Anshul Kundaje. See "https://sites.google.com/site/anshulkundaje/projects/epigenomeroadmap#TOC-Uniformly-processed-peak-calls" for more in-depth explanation of this pipeline.
# NOTE: Compares a single chip-seq dataset against a single control.
# =====================================================================
cell=$1

# Prefixes:
CPREF=FINAL_${EPITOPE}_${cell}
QCPREFIX=${QCDIR}/${CPREF}
OFPREFIX=${TADIR}/${CPREF}
FINAL_SUB_FILE="${OFPREFIX}.sub.tagAlign.gz" # Already subsampled through SCCA
FRAGLENFILE=${QCPREFIX}.sub.cc.qc

# Arguments / Software:
GENOMESIZE='hs' 
# Check if MACS2.1 is in path: 
# (load with reuse .macs2-2.1.1.20160309-python-2.7.1-sqlite3-rtrees)
if [[ -z $(which macs2) ]]; then echo 'ERROR: MACS executable not in $PATH'; exit 1; fi
# Alternative MACS2.0 path: MACSPATH='/broad/compbio/anshul/projects/encode/preprocessing/peakcalling/macs/scripts/macs2'

# ========================
# Data and Pre-processing:
# ========================
# Experiment data:
chipstub="${CPREF}.sub"
chip=${FINAL_SUB_FILE}

# Control data:
WCEDIR="${CHPDIR}/files/WCE"
controlstub="FINAL_WCE_${cell}"
control="${WCEDIR}/tagAlign/${controlstub}.tagAlign.gz"

# Provide for no control:
noControl=0
if [[ ! -s ${control} ]]
then
    noControl=1
    controlstub="noControl"
fi

# Output file name
PKPREF="${PKDIR}/${chipstub}_VS_${controlstub}"
if [[ -s "${PKPREF}_peaks.narrowPeak" ]]
then
    echo "DONE: Skipping ${PKPREF} as it is already done."
    exit 1
fi

# Get fraglen corresponding to ChIP file
if [[ -s ${FRAGLENFILE} ]]
then
    fraglen=$( awk '{printf("%0.0d",$3 / 2)}' ${FRAGLENFILE} )
    echo "Fragment length = $fraglen * 2"
else
    echo "WARNING: No fragment length found corresponding to file ${chip}"
fi

# ==========================
# Get read depth ratio and
# create tmp copy of control
# ==========================
tmpdir="${TMP}/tmp${RANDOM}_${RANDOM}"
mkdir ${tmpdir}

# Tmp unzipped version for macs:
combchip="${tmpdir}/${chipstub}.tagAlign"
if [[ -f "${combchip}" ]]; then rm -rf "${combchip}"; fi
gunzip -c ${chip} >> ${combchip}

echo "Counting ChIP reads: ${chip}"
chipReads=$(wc -l ${combchip} | awk '{printf "%f", $1/1000000}') 
sval=${chipReads}

if [[ ${noControl} == '0' ]]
then
    combcontrol="${tmpdir}/${controlstub}.tagAlign"
    if [[ -f "${combcontrol}" ]]; then rm -rf "${combcontrol}"; fi

    echo "Combining Control replicates: ${control}"
    gunzip -c ${control} >> ${combcontrol}
    controlReads=$(wc -l ${combcontrol} | awk '{printf "%f", $1/1000000}' )
    sval=$(echo "${chipReads} ${controlReads}" | awk '$1>$2{printf "%f",$2} $1<=$2{printf "%f",$1}' )
fi

# =========================
# NarrowPeak Calling (MACS)
# =========================
if [[ ${noControl} == '0' ]]
then
    if [[ ${FRAGLENFILE} == 'NULL' || ${fraglen} -eq 0 ]]
    then
        macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF} -g ${GENOMESIZE} -p 1e-2 --nomodel --extsize 73
    else
        macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF} -g ${GENOMESIZE} -p 1e-2 --nomodel --extsize ${fraglen}
    fi
else
    echo "ERROR: No Control provided; cannot run MACS2.1"
    exit 1
fi

rm -f ${PKPREF}_peaks.xls ${PKPREF}_peaks.bed ${PKPREF}_summits.bed

# Add rank and compress:
${BINDIR}/addRank2NarrowPeakFile.sh ${PKPREF}_peaks.narrowPeak 8 | gzip -c > ${PKPREF}.narrowPeak.gz
rm -f ${PKPREF}_peaks.narrowPeak

# ============================
# BroadPeak/GappedPeak Calling
# ============================
if [[ ${noControl} == '0' ]]
then
    if [[ ${FRAGLENFILE} == 'NULL' || ${fraglen} -eq 0 ]]
    then
        macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF}.broad -g ${GENOMESIZE} -p 1e-2 --broad --nomodel --extsize 73
    else
        macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF}.broad -g ${GENOMESIZE} -p 1e-2 --broad --nomodel --extsize ${fraglen}
    fi
else
    echo "ERROR: No Control provided; cannot run MACS2.1"
    exit 1
fi

# Compress:
gzip -c ${PKPREF}.broad_peaks.broadPeak > ${PKPREF}_broad.broadPeak.gz
gzip -c ${PKPREF}.broad_peaks.gappedPeak > ${PKPREF}_broad.gappedPeak.gz

# Clean up:
rm -rf ${tmpdir}
rm -f ${PKPREF}.broad_peaks.xls ${PKPREF}.broad_summits.bed ${PKPREF}.broad_peaks.broadPeak ${PKPREF}.broad_peaks.gappedPeak

