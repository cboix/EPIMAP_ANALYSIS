#!/bin/bash -l
# =====================================================================
# Call peaks and generate signal for single replicates
# (use this for resampling pipeline)
#   Adapted from Anshul Kundaje. See "https://sites.google.com/site/anshulkundaje/projects/epigenomeroadmap#TOC-Uniformly-processed-peak-calls" for more in-depth explanation of this pipeline.
# =====================================================================
chip=$1
chipstub=$2
controlstub=$3
QCPREFIX=$4
OUTDIR=$5
TYPE=$6

# Prefixes:
CPREF=FINAL_${EPITOPE}_${cell}
FRAGLENFILE=${QCPREFIX}.cc.qc

# Arguments / Software:
GENOMESIZE='hs' 
# Check if MACS2.1 is in path: 
# (load with reuse .macs2-2.1.1.20160309-python-2.7.1-sqlite3-rtrees)
if [[ -z $(which macs2) ]]; then echo 'ERROR: MACS executable not in $PATH'; exit 1; fi
# Alternative MACS2.0 path: MACSPATH='/broad/compbio/anshul/projects/encode/preprocessing/peakcalling/macs/scripts/macs2'

# ========================
# Data and Pre-processing:
# ========================
# Control data:
WCEDIR="${CHPDIR}/files/WCE"
control="${WCEDIR}/tagAlign/${controlstub}.tagAlign.gz"

# Output file name
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

# Tmp unzipped versions for macs:
combchip="${tmpdir}/${chipstub}.tagAlign"
if [[ -f "${combchip}" ]]; then rm -rf "${combchip}"; fi
gunzip -c ${chip} >> ${combchip}

combcontrol="${tmpdir}/${controlstub}.tagAlign"
if [[ -f "${combcontrol}" ]]; then rm -rf "${combcontrol}"; fi
echo "Combining Control replicates: ${control}"
gunzip -c ${control} >> ${combcontrol}

echo "Counting ChIP reads:"
chipReads=$(wc -l ${combchip} | awk '{printf "%f", $1/1000000}') 
sval=${chipReads}
controlReads=$(wc -l ${combcontrol} | awk '{printf "%f", $1/1000000}' )
sval=$(echo "${chipReads} ${controlReads}" | awk '$1>$2{printf "%f",$2} $1<=$2{printf "%f",$1}' )

# ===================
# Peak Calling (MACS)
# ===================
PKPREF="${OUTDIR}/${chipstub}_VS_${controlstub}"
if [[ "$TYPE" == "narrow" ]] 
then
    # Narrow:
    if [[ ${FRAGLENFILE} == 'NULL' || ${fraglen} -eq 0 ]]
    then
        macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF} -g ${GENOMESIZE} -p 1e-2 --nomodel --extsize 73
    else
        macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF} -g ${GENOMESIZE} -p 1e-2 --nomodel --extsize ${fraglen}
    fi

    rm -f ${PKPREF}_peaks.xls ${PKPREF}_peaks.bed ${PKPREF}_summits.bed

    # Add rank and compress:
    ${BINDIR}/addRank2NarrowPeakFile.sh ${PKPREF}_peaks.narrowPeak 8 | gzip -c > ${PKPREF}.narrowPeak.gz
    rm -f ${PKPREF}_peaks.narrowPeak

else 
    # BroadPeak/GappedPeak
    if [[ ${FRAGLENFILE} == 'NULL' || ${fraglen} -eq 0 ]]
    then
        macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF}.broad -g ${GENOMESIZE} -p 1e-2 --broad --nomodel --extsize 73
    else
        macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF}.broad -g ${GENOMESIZE} -p 1e-2 --broad --nomodel --extsize ${fraglen}
    fi

    # Compress:
    gzip -c ${PKPREF}.broad_peaks.broadPeak > ${PKPREF}_broad.broadPeak.gz
    gzip -c ${PKPREF}.broad_peaks.gappedPeak > ${PKPREF}_broad.gappedPeak.gz

    # Clean up:
    rm -f ${PKPREF}.broad_peaks.xls ${PKPREF}.broad_summits.bed ${PKPREF}.broad_peaks.broadPeak ${PKPREF}.broad_peaks.gappedPeak
fi

# ==================
# Signal generation:
# ==================
PKPREF="${OUTDIR}/${chipstub}_VS_${controlstub}_bdg"
if [[ ${FRAGLENFILE} == 'NULL' || ${fraglen} -eq 0 ]]
then
    macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF} -g ${GENOMESIZE} -p 1e-2 --nomodel --extsize 73 -B --SPMR
else
    macs2 callpeak -t ${combchip} -c ${combcontrol} -f BED -n ${PKPREF} -g ${GENOMESIZE} -p 1e-2 --nomodel --extsize ${fraglen} -B --SPMR
fi

# Remove unecessary files and compress files:
rm -rf ${tmpdir}
rm -f ${PKPREF}_peaks.xls ${PKPREF}_peaks.bed ${PKPREF}_summits.bed
gzip ${PKPREF}_peaks.narrowPeak

# ===================
# Generate bedgraphs:
# ===================
# NOTE: We do not generate the foldchange bedgraph here!
# -log10pval (only this is needed for dist. comparison)
LVPREF=${PKPREF}.pval.signal
if [[ ! -e "${LVPREF}.bigwig" ]]
then
    echo "Generating log10pval bedgraph"
    macs2 bdgcmp -t ${PKPREF}_treat_pileup.bdg -c ${PKPREF}_control_lambda.bdg -o ${PKPREF}_ppois.bdg -m ppois -S ${sval}
    slopBed -i ${PKPREF}_ppois.bdg -g ${CHROMSIZES} -b 0 | bedClip stdin ${CHROMSIZES} ${LVPREF}.bedgraph
    rm -rf ${PKPREF}_ppois.bdg

    echo "Sort and create log10pval bigWig file"
    sort -k1,1 -k2,2n ${LVPREF}.bedgraph > ${LVPREF}.bedgraph.tmp
    bedGraphToBigWig ${LVPREF}.bedgraph.tmp ${CHROMSIZES} ${LVPREF}.bigwig
    gzip -c ${LVPREF}.bedgraph.tmp > ${LVPREF}.bedgraph.gz
    rm -f ${LVPREF}.bedgraph.tmp ${LVPREF}.bedgraph
fi
rm -f ${PKPREF}_treat_pileup.bdg ${PKPREF}_control_lambda.bdg

