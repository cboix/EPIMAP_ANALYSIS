#!/bin/bash 
# ===============================================
# Script to run single files through the pipeline
# ===============================================
export EPITOPE=single_files
export FILEINFO=${CHPLNK}/${EPITOPE}_bam.tsv
export EPDIR=${CHPDIR}/files/${EPITOPE}
export TADIR=${EPDIR}/tagAlign
export PKDIR=${EPDIR}/peaks
export BDGDIR=${EPDIR}/bedgraph
export QCDIR=${EPDIR}/qc
mkdir -p $EPDIR $TADIR $QCDIR $BDGDIR $PKDIR ${CONVERTED_DATADIR}/${EPITOPE}
echo "Processing epitope: $EPITOPE - has $( wc -l $FILEINFO | awk '{print $1 - 1}' ) files"

# For example: 
SINGLEINFO=${ANNDIR}/single_info.tsv
echo "GTEX_N7MT_0011_R10a\t/broad/hptmp/zhizhuo/egtex_chipseq/nick24pilot_GTEX_N7MT_0011_R10a/combine4lane.bam\tzhizhuo\thg19" > $SINGLEINFO
echo "GTEX_R55F_0011_R10a\t/broad/hptmp/zhizhuo/egtex_chipseq/nick24pilot_GTEX_R55F_0011_R10a/combine4lane.bam\tzz2\thg19" >> $SINGLEINFO

# ====================================
# STEP1 - Download + Process each one:
# ====================================
# link rep assembly processing cell file
# Get attributes:
while read id location cell assembly 
do 
    CELL=$cell
    processing="unfiltered_alignments"
    # Pre-process to tagAlign:
    STEP1_FILE="${TADIR}/FINAL_${id}_${cell}.tagAlign.gz"
    if LC_ALL=C gzip -l ${STEP1_FILE} | awk 'NR==2 {exit($2!=0)}'
    then
        echo "- STEP1 Preprocess $id $cell"
        echo "- Making ${STEP1_FILE}"
        source $BINDIR/BAMtoTA_download_preprocess.sh $id $cell $location $assembly $processing
    fi

    # ==========================================
    # STEP2 - Pool same cell type + perform SCCA
    # Output a bedgraph for subsampled + not
    # ==========================================
    STEP2_FILE="${TADIR}/FINAL_${EPITOPE}_${CELL}.sub.tagAlign.gz"
    NREADS=30000000 # subsampling depth

    # Prefixes:
    CPREF=FINAL_${EPITOPE}_${cell}
    QCPREFIX=${QCDIR}/${CPREF}
    OFPREFIX=${TADIR}/${CPREF}
    FINAL_TA_FILE="${OFPREFIX}.tagAlign.gz"
    FINAL_SUB_FILE="${OFPREFIX}.sub.tagAlign.gz"

    # NOTE: all this also in Pool_SCCA.
    if LC_ALL=C gzip -l ${FINAL_TA_FILE} | awk 'NR==2 {exit($2!=0)}'
    then
        echo "Pool files"
        zcat ${STEP1_FILE} | gzip -c > ${FINAL_TA_FILE}

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
    if [[ ! -s ${QCPREFIX}.sub.cc.qc ]] 
    then
        # On subset only (for now)
        Rscript ${SPPDIR}/run_spp_nodups.R -c=${FINAL_SUB_FILE} -filtchr=chrM -savp=${QCPREFIX}.sub.cc.plot.pdf -out=${QCPREFIX}.sub.cc.qc
        sed -i -r 's/,[^\t]+//g' ${QCPREFIX}.sub.cc.qc

        echo "QC stats on the final tagAlign files."
        Rscript ${SPPDIR}/run_spp_nodups.R -rf -c=${FINAL_SUB_FILE} -savp=${QCPREFIX}.sub.tagAlign.pdf -out=${QCPREFIX}.sub.stats.qc
    fi

    # ===================================================
    # STEP3 - Peak + broadPeak calling: (WITH NO CONTROL)
    # ===================================================
    chipstub="FINAL_${EPITOPE}_${CELL}.sub"
    CPREF=${chipstub}
    outPref="${PKDIR}/${CPREF}" 
    STEP3_FILE="${outPref}_peaks.narrowPeak.gz"
    if [[ "$( LC_ALL=C gzip -l ${STEP3_FILE} | awk 'NR==2 {print $2}' )" == "0" ]] 
    then
        echo "- STEP3 Call narrow and broad peaks for ${CELL} in ${EPITOPE}"
        echo "- Making $STEP3_FILE"
        chip=${FINAL_TA_FILE}

        # Prefixes:
        CPREF=FINAL_${EPITOPE}_${cell}
        FRAGLENFILE=${QCPREFIX}.sub.cc.qc
        GENOMESIZE='hs' 

        fraglen=$( awk '{printf("%0.0d",$3 / 2)}' ${FRAGLENFILE} )
        echo "Fragment length = $fraglen * 2"

        # Create tmp copy of chip
        tmpdir="${TMP}/tmp${RANDOM}_${RANDOM}"
        mkdir ${tmpdir}

        # Tmp unzipped versions for macs:
        combchip="${tmpdir}/${chipstub}.tagAlign"
        if [[ -f "${combchip}" ]]; then rm -rf "${combchip}"; fi
        gunzip -c ${chip} >> ${combchip}

        # Peak Calling (MACS)
        PKPREF="${PKDIR}/${chipstub}"
        # Narrow:
        if [[ ${FRAGLENFILE} == 'NULL' || ${fraglen} -eq 0 ]]
        then
            macs2 callpeak -t ${combchip} -f BED -n ${PKPREF} -g ${GENOMESIZE} -p 1e-2 --nomodel --extsize 73
        else
            macs2 callpeak -t ${combchip} -f BED -n ${PKPREF} -g ${GENOMESIZE} -p 1e-2 --nomodel --extsize ${fraglen}
        fi

        rm -f ${PKPREF}_peaks.xls ${PKPREF}_peaks.bed ${PKPREF}_summits.bed

        # Add rank and compress:
        ${BINDIR}/addRank2NarrowPeakFile.sh ${PKPREF}_peaks.narrowPeak 8 | gzip -c > ${PKPREF}.narrowPeak.gz
        rm -f ${PKPREF}_peaks.narrowPeak
    fi

done < $SINGLEINFO
