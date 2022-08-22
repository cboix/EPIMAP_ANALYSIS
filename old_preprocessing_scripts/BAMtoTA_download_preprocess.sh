#!/bin/bash
# Pipeline for individual DNase-seq runs:
SAMTL="/broad/software/free/Linux/redhat_6_x86_64/pkgs/samtools/samtools_1.3/bin/samtools" # Requires the 1.3 version.

id=$1;
cell=$2;
link=$3;
assembly=$4;
processing=$5;

RLEN=36

# Skip quality filter if file is filtered.
if [[ "$processing" == "alignments" ]] 
then
    echo "Setting no MAPQ threshold for alignments because they have already been filtered. Some filtered files are Phred64, which causes problems with using the consistent MAPQ=30."
    MAPQ_THRESH=0
else
    echo "Setting MAPQ threshold at 30."
    MAPQ_THRESH=30
fi


# Prefixes:
CPREF="${id}_${cell}_${assembly}"
PREFIX=${TADIR}/${CPREF}
QC_PREFIX=${QCDIR}/${CPREF}

# Files
RAW_BAM_FILE=${PREFIX}.bam
SORTED_BAM=${PREFIX}.sorted.bam
QUAL_QC="${QC_PREFIX}.mapq.qc" # MAPQ range in raw file

# Filtering files: 
FILT_BAM_PREFIX="${PREFIX}.filt.srt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
DUP_FILE_QC="${QC_PREFIX}.filt.srt.dup.qc" # QC file for marking duplicates

# Final BAM files:
FINAL_QC_PREFIX="${QC_PREFIX}.filt.nodup.srt"
FINAL_BAM_PREFIX="${PREFIX}.filt.nodup.srt"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" # To be stored
FINAL_BAM_FILE_MAPSTATS="${FINAL_QC_PREFIX}.flagstat.qc" # QC file
PBC_FILE_QC="${FINAL_QC_PREFIX}.pbc.qc" # quality stats

# Final TA files:
FINAL_TA_FILE="${FINAL_BAM_PREFIX}.SE.tagAlign.gz"
OFNAME=${TADIR}/FINAL_${id}_${cell}.tagAlign

# PIPELINE:
if LC_ALL=C gzip -l ${OFNAME} | awk 'NR==2 {exit($2!=0)}'
then
    echo "Empty/Missing final TagAlign file"

    LINKEXT=${link#*.*.*.}
    FULLEXT=${link##*.}
    echo "Link extension is of the - $LINKEXT - filetype."
    if [[ "$LINKEXT" == "tagAlign.gz" ]] 
    then 
        # STARTING FROM TAGALIGN:
        echo "Downloading ${FINAL_TA_FILE}"
        curl -L -u 'R7GOMSZL:7sccfx4rtdb6sqse' -o ${FINAL_TA_FILE} ${link} # submitter key
    elif [[ "$FULLEXT" == "bam" ]] 
    then
        # =======================
        # STARTING FROM BAM FILES
        # =======================
        if $SAMTL view -F 0x904 -c ${RAW_BAM_FILE} | awk '{exit($1!=0)}'
        then 
            if [[ "$(ls -s $link | awk '$1 > 0{print 1}' )" == "1" ]] 
            then
                cp $link ${RAW_BAM_FILE}
                # $SAMTL sort $link > ${RAW_BAM_FILE}
                echo "Indexing ${RAW_BAM_FILE}"
                $SAMTL index ${RAW_BAM_FILE} # Index
            else
                # Download data:
                if $SAMTL view -F 0x904 -c ${RAW_BAM_FILE} | awk '{exit($1!=0)}'
                then 
                    echo "Downloading ${RAW_BAM_FILE}"
                    curl -L -u 'R7GOMSZL:7sccfx4rtdb6sqse' -o ${RAW_BAM_FILE} ${link} # submitter key
                    echo "Indexing ${RAW_BAM_FILE}"
                    $SAMTL index ${RAW_BAM_FILE} # Index
                fi
            fi
        fi
        # TODO FIGURE OUT WHAT VERSION OF ILLUMINA EACH HAS
        # Store first 500 lines for diagnostics:
        $SAMTL view ${RAW_BAM_FILE} | head -n 500 > ${PREFIX}_f500.sample.sam
        # TODO What to do with 50/100 bp reads etc? 

        #================================================
        # Filter file by read quality and mark duplicates 
        #================================================
        if $SAMTL view -F 0x904 -c ${FILT_BAM_FILE} | awk '{exit($1!=0)}'
        then
            # =========================================================================================
            # Remove with 1804: unmapped, mate unmapped, not primary alignment, reads failing platform
            # Remove low MAPQ reads
            # Take a look at SAM flags here: https://broadinstitute.github.io/picard/explain-flags.html
            # =========================================================================================
            # MAPQ quality scores: aka do we filter by MAPQ 30? 
            if [[ ! -s ${QUAL_QC} ]] 
            then
                $SAMTL view ${RAW_BAM_FILE} | cut -f5 | sort -n | uniq -c > ${QUAL_QC}
            fi

            echo "Filtering to get ${FILT_BAM_FILE}"
            # Ensure that $SAMTL v1.3 is used.
            $SAMTL view -F 1804 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} | $SAMTL sort - -T ${TMP}/${CPREF}.tmp -o ${FILT_BAM_FILE}
            $SAMTL view -H ${FILT_BAM_FILE} | grep SO

            echo "Marking duplicates for ${FILT_BAM_FILE}"
            TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
            MARKDUP="/seq/software/picard/1.802/bin/MarkDuplicates.jar";
            #MARKDUP="/seq/software/picard/current/bin/MarkDuplicates.jar"; # This breaks.

            java -Xmx4G -jar ${MARKDUP} INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false # Crashes if low memory

            mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}
        fi

        # =============================================
        # Remove dups & Index final position sorted BAM
        # =============================================
        if $SAMTL view -F 0x904 -c ${FINAL_BAM_FILE} | awk '{exit($1!=0)}'
        then
            echo "Removing duplicates to get ${FINAL_BAM_FILE}"
            $SAMTL view -F 1804 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE} # NOTE 1804 is for unpaired reads!

            echo "Indexing ${FINAL_BAM_FILE}"
            $SAMTL index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}

            $SAMTL flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}

            # =============================
            # Compute library complexity
            # =============================
            # Obtain unique count statistics
            # PBC File output:
            # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
            echo "== COMPUTE LIBRARY COMPLEXITY =="
            bedtools bamtobed -i ${FILT_BAM_FILE} | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | \
                awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\n",mt,m0,m1,m2}' > ${PBC_FILE_QC}
        fi

        # Remove temporary data if previous steps have non-zero output
        if $SAMTL view -F 0x904 -c ${FINAL_BAM_FILE} | awk '{exit($1==0)}'
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
            # cut all artificial reads to RLEN (36 - match mappability track).
            echo "== Make TAGALIGN file =="
            chrpref=$( samtools view ${FINAL_BAM_FILE} | head | awk 'NR==1{sub("[0-9]+","",$3); out=($3=="chr"?"":"chr"); print out}' ) # ensure chrSpace for nums
            bedtools bamtobed -i ${FINAL_BAM_FILE} | awk -v chrpref=$chrpref 'BEGIN{OFS="\t"}{$1=chrpref$1; $4="N";$5="1000";print $0}' | \
                awk 'BEGIN{FS="\t";OFS="\t"} $6=="+"{$3=$2+'"${RLEN}"';print $0} $6=="-"{$2=$3-'"${RLEN}"';print $0}' | \
                gzip -c > ${FINAL_TA_FILE}
        fi
    else 
        echo "Unrecognized filetype: $LINKEXT or $FULLEXT"
        break
    fi

    # ============================
    # LiftOver to hg38 before UMAP
    # ============================
    if [[ "$assembly" == "hg19" ]] 
    then
        echo "== LIFTOVER hg19 to hg38 =="
        TMPOVER=${PREFIX}.over.tagAlign
        liftOver ${FINAL_TA_FILE} ${LOCHAIN} ${TMPOVER} ${PREFIX}.over.unmapped
        gzip -c ${TMPOVER} > ${FINAL_TA_FILE}
        rm ${TMPOVER}
    fi

    ###################################################
    ### mapFilterTagAlignFiles.sh
    ### THIS CODE IS BASED ON ANSHUL'S mapFilterTagAlignFiles.sh CODE.
    ### IT CALLS filterUniqueReads TO REMOVE ANY POSSIBLE READS IN 'UNMAPPABLE' LOCATIONS
    ###################################################
    # REMOVE UNMAPPABLE REGIONS:
    if LC_ALL=C gzip -l ${OFNAME} | awk 'NR==2 {exit($2!=0)}'
    then
        echo "== Remove hg38 UMAP regions =="
        bedtools intersect -a ${FINAL_TA_FILE} -b ${MAP36} | gzip -c > ${OFNAME}.gz
        # Sanity check to see if we were successful:
        zcat ${OFNAME} | tail
    fi

    if LC_ALL=C gzip -l ${OFNAME} | awk 'NR==2 {exit($2==0)}'
    then
        rm -rf ${FINAL_BAM_FILE}
        rm -rf ${FINAL_TA_FILE}
    fi
fi

# Print out number of reads: 
zcat ${OFNAME}.gz | wc -l > ${QC_PREFIX}.numreads

