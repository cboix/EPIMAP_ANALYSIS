#!/bin/bash
# Pipeline for individual ATAC-seq runs:

# TODO READ IN INFORMATION 
eval $( sed -n "$SGE_TASK_ID"p $INFO | awk 'BEGIN{FS="\t"}{printf("EXPID=%s; EXPNUM=%s; CONDITION=\"%s\"; CELL=%s;",$1,$2,$3,$4)}' )

MAPQ_THRESH=30

# PREFIXES AND FILES: 
PREFIX=$BAMDIR/${DATAPREF}_${EXPNUM}

RAW_BAM_FILE=${PREFIX}_mm9_bestmap.bam
SORTED_BAM=${PREFIX}_mm9_bestmap.sorted.bam
# QUALBAM=${PREFIX}.best.qual.bam

# Filtering files: 
FILT_BAM_PREFIX="${PREFIX}.filt.srt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc" # QC file for marking duplicates

# Final BAM files:
FINAL_BAM_PREFIX="${PREFIX}.filt.nodup.srt"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" # To be stored
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file
PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc" # quality stats

# Final TA files:
FINAL_TA_FILE="${FINAL_BAM_PREFIX}.SE.tagAlign.gz"
FINAL_TA_MAP_FILE="${FINAL_BAM_PREFIX}.SE.map.tagAlign.gz"

# Set output and log file name (FIXME not necessarily best place to define these)
OFNAME=$(echo ${FINAL_TA_FILE} | sed -r -e 's/\.tagAlign\.gz$/.map.tagAlign/g')
LOGFILE=$(echo ${FINAL_TA_FILE} | sed -r -e 's/\.tagAlign\.gz$/.map.tagAlign.logfile/g')

# PIPELINE:
if LC_ALL=C gzip -l ${OFNAME} | awk 'NR==2 {exit($2!=0)}'
then

    #================================================
    # Filter file by read quality and mark duplicates 
    #================================================
    # TODO add filter for mapped tagAlign
    if samtools view -F 0x904 -c ${FILT_BAM_FILE} | awk '{exit($1!=0)}'
    then
        # =========================================================================================
        # Remove with 1804: unmapped, mate unmapped, not primary alignment, reads failing platform
        # Remove low MAPQ reads
        # Take a look at SAM flags here: https://broadinstitute.github.io/picard/explain-flags.html
        # =========================================================================================
        echo "Filtering to get ${FILT_BAM_FILE}"
        samtools view -F 1804 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} | samtools sort - -T ${FILT_BAM_PREFIX}.tmp -o ${FILT_BAM_FILE}
        samtools view -H ${FILT_BAM_FILE} | grep SO

        echo "Marking duplicates for ${FILT_BAM_FILE}"
        TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
        MARKDUP="/seq/software/picard/1.802/bin/MarkDuplicates.jar";
        #MARKDUP="/seq/software/picard/current/bin/MarkDuplicates.jar"; # might break

        java -Xmx4G -jar ${MARKDUP} INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false # Crashes if low memory

        mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}
    fi

    # =============================================
    # Remove dups & Index final position sorted BAM
    # =============================================
    if samtools view -F 0x904 -c ${FINAL_BAM_FILE} | awk '{exit($1!=0)}'
    then
        echo "Removing duplicates to get ${FINAL_BAM_FILE}"
        samtools view -F 1804 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE} # NOTE 1804 is for unpaired reads!

        echo "Indexing ${FINAL_BAM_FILE}"
        samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}

        samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}

        # =============================
        # Compute library complexity
        # =============================
        # Obtain unique count statistics
        # PBC File output:
        # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
        bedtools bamtobed -i ${FILT_BAM_FILE} | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | \
            awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\n",mt,m0,m1,m2}' > ${PBC_FILE_QC}
    fi

    # Remove temporary data if previous steps have non-zero output
    if samtools view -F 0x904 -c ${FINAL_BAM_FILE} | awk '{exit($1==0)}'
    then
        rm ${FILT_BAM_FILE}
        rm ${RAW_BAM_FILE} ${RAW_BAM_FILE}.bai
    fi

    ###################################################
    ### STEP 2A --- BASED ON ANSHUL'S ENCODE3 PROPOSAL
    ###################################################

    # ===================
    # Create tagAlign file
    # ===================
    if LC_ALL=C gzip -l ${FINAL_TA_FILE} | awk 'NR==2 {exit($2!=0)}'
    then
        bedtools bamtobed -i ${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | \
            awk 'BEGIN{FS="\t";OFS="\t"} $6=="+"{$3=$2+'"${RLEN}"';print $0} $6=="-"{$2=$3-'"${RLEN}"';print $0}' | \
            gzip -c > ${FINAL_TA_FILE}
    fi

    ###################################################
    ### mapFilterTagAlignFiles.sh
    ### THIS CODE IS BASED ON ANSHUL'S mapFilterTagAlignFiles.sh CODE.
    ### IT CALLS filterUniqueReads TO REMOVE ANY POSSIBLE READS IN 'UNMAPPABLE' LOCATIONS
    ###################################################

    # TODO having problems with MATLAB!
    if LC_ALL=C gzip -l ${OFNAME} | awk 'NR==2 {exit($2!=0)}'
    then
        export TMPLOC=${TMP}/tmp_${RANDOM}${RANDOM};
        mkdir -p ${TMPLOC};
        export MCR_CACHE_ROOT=${TMPLOC};
        ${FILTDIR}/filterUniqueReads -s=${SDIR} -u=${UDIR} -v=${LOGFILE} ${FINAL_TA_FILE} | grep -v 'Warning' | gzip -c > ${OFNAME}.gz
        rm -rf ${TMPLOC}

        # Sanity check to see if we were successful:
        zcat ${OFNAME} | tail
    fi

    # INSTEAD RUN: 
    if [[ "0" == "1" ]]
    then

        export TMPLOC=${TMPDIR}/tmp_${RANDOM}${RANDOM};
        mkdir -p ${TMPLOC};
        export MCR_CACHE_ROOT=${TMPLOC};

        while read TAFILE
        do
            OFNAME=$(echo ${TAFILE} | sed -r -e 's/\.tagAlign\.gz$/.map.tagAlign/g')
            LOGFILE=$(echo ${TAFILE} | sed -r -e 's/\.tagAlign\.gz$/.map.tagAlign.logfile/g')

            echo "$TAFILE to $OFNAME\.gz"
            ${FILTDIR}/filterUniqueReads -s=${SDIR} -u=${UDIR} -v=${LOGFILE} ${TAFILE} | grep -v 'Warning' | gzip -c > ${OFNAME}.gz


            # Sanity check to see if we were successful:
            zcat ${OFNAME} | tail
        done < <(ls $BAMDIR/*SE.tagAlign.gz)

        rm -rf ${TMPLOC}

    fi


    if LC_ALL=C gzip -l ${OFNAME} | awk 'NR==2 {exit($2==0)}'
    then
        rm ${FINAL_BAM_FILE}
        rm ${FINAL_TA_FILE}
    fi

fi

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Peak calling and visualization:
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create wigfiles to view data:

BWPREF=${BWDIR}/${EXPID}
PKPREF=${PKDIR}/${EXPID}_peaks

# FIXME into separate directory:
if [[ ! -s ${BWPREF}.only.bw ]] 
then
    # Unzip, sort, and then convert to a bedgraph file
    zcat ${OFNAME}.gz | sort -k1,1 -k2,2n | genomeCoverageBed -bg -i stdin -g ${CHROMSIZES} > ${BWPREF}.bedGraph

    # Convert the bed graph to a bigWig file for viewing 
    bedGraphToBigWig ${BWPREF}.bedGraph ${CHROMSIZES} ${BWPREF}.only.bw
fi

PKPREF=${PKDIR}/${EXPID}
cd ${PKDIR}
if [[ ! -s ${PKPREF}_peaks.narrowPeak ]]
then
    # $SFTDIR/macs2 callpeak -t ${FINAL_BAM_FILE} -f BAMPE -n ${EXPID} --nomodel  --shift -100 --extsize 200 --gsize 1.87e9
    $SFTDIR/macs2 callpeak -t ${OFNAME}.gz -f BED -n ${EXPID} --nomodel  --shift -100 --extsize 200 --gsize 1.87e9

    # Filter MACS peaks by q-value (something like the following):
    awk -v eid=$EXPID 'BEGIN{OFS="\t"; FS="\t"}$9 >= 5{curRow=curRow+1; print $0,eid curRow}' ${PKPREF}_peaks.narrowPeak > ${PKPREF}_peaks.filt.bed 
fi

PKBIN=${PKPREF}_filt.gz
PKINT=${PKPREF}_filt_logQ.gz
if [[ ! -s $PKBIN || ! -s $PKINT ]]
then
    bedtools sort -i ${PKPREF}_peaks.filt.bed > ${PKPREF}.tmp.bed

    # Intersects into binned track:
    bedtools intersect -a $WINDOWS -b ${PKPREF}.tmp.bed -loj -wao -sorted | awk '{a=($12>0?$12:0); print $1,$2,$3,a}' > ${PKPREF}_over.tmp

    awk -v expid=$EXPID 'BEGIN{OFS="\t"; f="chr\tstart\tend"; a = expid}{e = sprintf("%s\t%s\t%s",$1,$2,$3); if(f==e){a=(a>$4?a:$4)} else {print a; a = $4}; f = e;}END{print a}' ${PKPREF}_over.tmp | gzip -c > ${PKINT}

    if [[ ! -s $PKBIN ]] 
    then
        bedtools intersect -a $WINDOWS -b ${PKPREF}.tmp.bed -loj -wa -c -sorted | awk -v expid=$EXPID 'BEGIN{print expid}{a=($4>0?1:0); print a}' | gzip -c > ${PKBIN}
    fi

    rm ${PKPREF}.tmp.bed
fi

ALLBIN=${PKDIR}/all_filt.tsv
FINALBIN=${PKDIR}/all_filt.reduced.tsv.gz
if [[ "0" == "1" ]]
then
    cp $WINDOWS $ALLBIN
    while read filtgz
    do
        echo $filtgz
        zcat $filtgz | pr -mts $ALLBIN - > ${ALLBIN}.tmp
        mv ${ALLBIN}.tmp ${ALLBIN}
    done < <(ls $PKDIR/*_filt.gz)

    awk '{a = 0; for(i=4;i<=NF;i++){a = a+ $i}; if(a > 0){print $0}}' ${ALLBIN} | gzip -c > ${FINALBIN}
    gzip $ALLBIN # finalize compression (1.2G to 80 M)!

fi

ALLINT=${PKDIR}/all_filt_logQ.tsv
FINALINT=${PKDIR}/all_filt_logQ.reduced.tsv.gz
if [[ "0" == "1" ]]
then
    cp $WINDOWS $ALLINT
    while read filtgz
    do
        echo $filtgz
        zcat $filtgz | pr -mts $ALLINT - > ${ALLINT}.tmp
        mv ${ALLINT}.tmp ${ALLINT}
    done < <(ls $PKDIR/*_filt_logQ.gz)

    awk '{a = 0; for(i=4;i<=NF;i++){a = a+ $i}; if(a > 0){print $0}}' ${ALLINT} | gzip -c > ${FINALINT}
    gzip $ALLINT # compress to interesting peaks
fi


# TODO make same as FINALINT but not filtered?


# while read prefix
# do
#     REG_PREF=${CEDIR}/${prefix}
#     BEDFILE=${REG_PREF}_hg19_peaks.bed.gz

#     if [[ ! -s ${REG_PREF}.nearest.gene.bed.gz ]] 
#     then
#         echo "Get TSS data for ${prefix}"
#         # Get closest transcript start site:
#         zcat $BEDFILE | bedtools sort -i - | bedtools closest -d -a - -b ${TSS_DATA} | gzip -c > ${REG_PREF}.tss.bed.gz
#         zcat ${REG_PREF}.tss.bed.gz | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$10,$11}' | sort -u | bedtools sort -i - | gzip -c > ${REG_PREF}.nearest.gene.bed.gz
#         # NOTE: If overlapping multiple genes, they are all the "closest"
#         # NOTE: All chr. M enhancers are -1 because there are no TSS in the gtf file in chrM

#     fi
# done < region_prefixes

# # TODO Look at distribution of proximity to TSS / call promoters.  (plotting of all?)

# # TODO merge replicates first.
# # Would like to cross-label enhancers in different tissues as the same? 

# # Make masterfile: 
# MASTER_PREF=$CEDIR/CE_MERGED
# MBED=${MASTER_PREF}_hg19_peaks.bed.gz

# # Assume MACS called indpt peaks
# # APPROACH 1: NAIVELY MERGE ANY OVERLAPS!
# # TODO BAD.
# while read prefix
# do
#     REG_PREF=${CEDIR}/${prefix}
#     BEDFILE=${REG_PREF}_hg19_peaks.bed.gz
#     if [[ -s $MBED ]] 
#     then
#         echo "Merging peaks from $prefix to master peaks file"
#         zcat $BEDFILE | bedtools sort -i - | bedtools closest -d -a - -b ${MBED} | awk '$11 == 0' > testing 

#         # Resolve all mappings and return a merged list!

#         > ${REG_PREF}.tss.bed.gz
#     else 
#         zcat $BEDFILE | awk '{sub("MACS","Master",$0); print $0}' | bedtools sort -i - | gzip -c > $MBED
#     fi

# done < region_prefixes




