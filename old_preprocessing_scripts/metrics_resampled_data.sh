#!/bin/bash
# SGE Array to process and generate metrics for resampled data.
SAMTL="/broad/software/free/Linux/redhat_6_x86_64/pkgs/samtools/samtools_1.3/bin/samtools" # Requires the 1.3 version.
SAMPQC=${SAMPDIR}/qc
SAMPPK=${SAMPDIR}/peaks
# ----------------
# Quality metrics (
# SCCA fraglength, NSC, and RSC
# PEAK calls - percent similarity
# Correlation of signal scores (using chromImpute globaldist function) 
# ----------------
export INFOLINE=$( sed "${SGE_TASK_ID}q;d" ${SAMPINFO} )
eval $( echo $INFOLINE | awk '{a=$1; sub(".tagAlign.gz","",a); sub(".*/","",a); printf("INFILE=%s; CELL=%s; TYPE=%s; PREF=%s",$1,$2,$3,a)}' )
QCPREFIX=${SAMPDIR}/qc/${PREF}

# ========================
# Compute quality metrics:
# ========================
if [[ ! -s ${QCPREFIX}.pbc.qc ]] 
then
    echo "Compute library complexity:"
    zcat ${INFILE} | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\n",mt,m0,m1,m2}' > ${QCPREFIX}.pbc.qc
fi

# ========================
# Use SPP for SCCA and QC:
# ========================
if [[ ! -s ${QCPREFIX}.cc.qc ]] 
then
    echo "Running SPP for SCCA on the ${CPREF} TA files."
    Rscript ${SPPDIR}/run_spp_nodups.R -c=${INFILE} -filtchr=chrM -savp=${QCPREFIX}.cc.plot.pdf -out=${QCPREFIX}.cc.qc
    sed -i -r 's/,[^\t]+//g' ${QCPREFIX}.cc.qc # Keeps best length.

    echo "QC stats on the final tagAlign files."
    Rscript ${SPPDIR}/run_spp_nodups.R -rf -c=${INFILE} -savp=${QCPREFIX}.tagAlign.pdf -out=${QCPREFIX}.stats.qc
fi

# ====================================
# Call peaks + generate signal tracks:
# ====================================
chipstub="${PREF}"
controlstub="FINAL_WCE_${CELL}"

CPREF=${chipstub}_VS_${controlstub}
outPref="${SAMPPK}/${CPREF}" 
STEP3_FILE="${outPref}_peaks.narrowPeak.gz"
if LC_ALL=C gzip -l ${STEP3_FILE} | awk 'NR==2 {exit($2!=0)}'
then
    echo "Call peaks and generate signal tracks:"
    # Shorter version of pipeline just for single rep peak calling: 
    source $BINDIR/single_rep_peaks_signal.sh ${INFILE} ${chipstub} ${controlstub} ${QCPREFIX} ${SAMPK} ${TYPE}

fi

# =============================================
# STEP4B - Convert pval files for ChromImpute: 
# =============================================
# STEP4B_FILE=${CONVERTED_DATADIR}/${EPITOPE}/chr20_${CPREF}.pval.signal.bedgraph.gz.wig.gz
# if LC_ALL=C gzip -l ${STEP4B_FILE} | awk 'NR==2 {exit($2!=0)}'
# then
#     echo "Convert signal into 25bp average tracks for ChromImpute"
#     echo "${cell}\t${EPITOPE}\t${chipstub}_VS_${controlstub}.pval.signal.bedgraph.gz" > ${outPref}_tab.txt
#     java -mx4000M -jar ${CHROMIMPUTE} Convert ${BDGDIR} ${PKPREF}_tab.txt ${CHROMSIZES} ${CONVERTED_DATADIR}/${EPITOPE}
#     rm ${PKPREF}_tab.txt
# fi
# TODO Implement IDR pipeline:
# ===========================
# STEP5 - Pseudo-reps and IDR
# ===========================
