#!/bin/bash -l
# =====================================================================
# Generate the signal track files in the ROADMAP/ENCODE processing pipeline (for imputation)
#   1. Fold change over control at base resolution
#   2. -log10p.val by Poisson model at base resolution
# 
#   Adapted from Anshul Kundaje. See "https://sites.google.com/site/anshulkundaje/projects/epigenomeroadmap#TOC-Uniformly-processed-and-normalized-genome-wide-signal-coverage-tracks" for more in-depth explanation of this pipeline.
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
PKPREF="${BDGDIR}/${chipstub}_VS_${controlstub}"
if [[ -s "${PKPREF}.pval.signal.bedgraph" ]]
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

# ==================================
# Script for signal track generation
# (Previously submitted separately).
# ==================================
tmpdir="${TMP}/tmp${RANDOM}_${RANDOM}"
mkdir ${tmpdir}

# --------------------------
# Get read depth ratio and
# create tmp copy of control
# --------------------------
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

# ===================================================================
# Create the pileup and control lambda bedgraph tracks using MACS2.1:
# ===================================================================
if [[ ${noControl} == '0' ]]
then
    if [[ ${FRAGLENFILE} == 'NULL' || ${fraglen} -eq 0 ]]
    then
        macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF} -g ${GENOMESIZE} -p 1e-2 --nomodel --extsize 73 -B --SPMR
    else
        macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF} -g ${GENOMESIZE} -p 1e-2 --nomodel --extsize ${fraglen} -B --SPMR
    fi
else
    # FIXME Provide an option for noControl runs!
    echo "ERROR: No Control provided; cannot run MACS2.1"
    exit 1
fi

# Remove unecessary files and compress files:
rm -rf ${tmpdir}
rm -f ${PKPREF}_peaks.xls ${PKPREF}_peaks.bed ${PKPREF}_summits.bed
gzip ${PKPREF}_peaks.narrowPeak

# ===================
# Generate bedgraphs:
# ===================
# FoldChange: 
FCPREF=${PKPREF}.fc.signal
if [[ ! -e "${FCPREF}.bedgraph.gz" ]]
then
    echo "Generating FoldChange bedgraph"
    macs2 bdgcmp -t ${PKPREF}_treat_pileup.bdg -c ${PKPREF}_control_lambda.bdg -o ${PKPREF}_FE.bdg -m FE
    slopBed -i ${PKPREF}_FE.bdg -g ${CHROMSIZES} -b 0 | bedClip stdin ${CHROMSIZES} ${FCPREF}.bedgraph
    rm -f ${PKPREF}_FE.bdg

    echo "Sort and create FC bigWig file"
    sort -k1,1 -k2,2n ${FCPREF}.bedgraph > ${FCPREF}.bedgraph.tmp
    # bedGraphToBigWig ${FCPREF}.bedgraph.tmp ${CHROMSIZES} ${LVPREF}.bigwig
    gzip -c ${FCPREF}.bedgraph.tmp > ${FCPREF}.bedgraph.gz
    rm -f ${FCPREF}.bedgraph.tmp ${FCPREF}.bedgraph
fi

# -log10pval: 
LVPREF=${PKPREF}.pval.signal
if [[ ! -e "${LVPREF}.bedgraph.gz" ]]
then
    echo "Generating log10pval bedgraph"
    macs2 bdgcmp -t ${PKPREF}_treat_pileup.bdg -c ${PKPREF}_control_lambda.bdg -o ${PKPREF}_ppois.bdg -m ppois -S ${sval}
    slopBed -i ${PKPREF}_ppois.bdg -g ${CHROMSIZES} -b 0 | bedClip stdin ${CHROMSIZES} ${LVPREF}.bedgraph
    rm -rf ${PKPREF}_ppois.bdg

    echo "Sort and create log10pval bigWig file"
    sort -k1,1 -k2,2n ${LVPREF}.bedgraph > ${LVPREF}.bedgraph.tmp
    # bedGraphToBigWig ${LVPREF}.bedgraph.tmp ${CHROMSIZES} ${LVPREF}.bigwig
    gzip -c ${LVPREF}.bedgraph.tmp > ${LVPREF}.bedgraph.gz
    rm -f ${LVPREF}.bedgraph.tmp ${LVPREF}.bedgraph
fi
rm -f ${PKPREF}_treat_pileup.bdg ${PKPREF}_control_lambda.bdg

