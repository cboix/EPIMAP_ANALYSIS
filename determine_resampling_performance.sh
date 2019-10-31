#!/bin/bash
# ============================================================
# Run downsampling and resampling on well profiled replicates 
# to assess resampling performance and potential biases.
# ============================================================
SAMPQC=${SAMPDIR}/qc
SAMPPK=${SAMPDIR}/peaks
mkdir -p $SAMPQC $SAMPPK

# ==========================================================
# List of tagAlign replicates with the following qualities: 
#   1. High read depth 
#   2. Range of SCCA quality scores
#   3. In Control, Histone ChIP, CR ChIP, and DNase
# ==========================================================
REPLICATES=1 # FIXME UP this number later
DOWN=( 5 10 20 ) 
UP=( 30 45 60 )

# Marks to use:
# -------------
# H3K27ac or H3K4me3 (for punctate)
# H3K27me3 or H3K36me3 (for broad)
export HIGHQDATA=${ANNDIR}/high_depth_resampling.tsv
export SAMPINFO=${SAMPDIR}/depth_resampling_files.tsv
rm -f ${SAMPINFO} # New lines will be added below:

IFS=$'\t'
while read CELL mark ID type assembly
do 
    REPLIST=$( seq 1 $REPLICATES )
    for REPL in ${REPLIST[@]}
    do 
        # Originals and prefixes:
        echo "$CELL $mark $type" 
        EPDIR=${CHPDIR}/files/${mark}
        TAFILE=${EPDIR}/tagAlign/FINAL_${ID}.tagAlign.gz
        QCPREF=${EPDIR}/qc/${ID}_${assembly} 
        SPREF=${SAMPDIR}/${ID}_${mark}_rep${REPL}

        # ----------------------------------------------
        # Generate files via subsampling and resampling:
        # ----------------------------------------------
        # Subsample: 
        for SUBS in ${DOWN[@]}
        do 
            echo "Subsampling to $SUBS"
            NREADS=${SUBS}000000 # depth in millions of reads
            PREFIX=${SPREF}_s${SUBS}
            SUBFILE=${PREFIX}.tagAlign.gz
            if [[ ! -s ${SUBFILE} ]] 
            then
                zcat ${TAFILE} | shuf -n ${NREADS} | sort -k1,1V -k2,2g | gzip -c > ${SUBFILE}
            fi
            echo "${SUBFILE} ${CELL} ${type}" >> $SAMPINFO

            # Resample up to higher depth:
            for RESAMP in ${UP[@]}
            do 
                if (( $SUBS < $RESAMP ))
                then
                    echo "Resampling to $RESAMP from $SUBS"
                    NREADS=${RESAMP}000000 
                    UPFILE=${PREFIX}_r${RESAMP}.tagAlign.gz
                    if [[ ! -s ${UPFILE} ]] 
                    then
                        zcat ${TAFILE} | shuf -n ${NREADS} -r | sort -k1,1V -k2,2g | gzip -c > ${UPFILE}
                    fi
                    echo "${UPFILE} ${CELL} ${type}" >> $SAMPINFO
                fi
            done
        done
    done

done < ${HIGHQDATA}

# Profile all files generated: (short queue + array jobs) 
SNUM=$(wc -l ${SAMPINFO} | awk '{print $1}')
qsub -cwd -t 1-$SNUM -l h_vmem=20G -N resamp_metrics -j y -b y -V -r y -o $DBDIR/out $BINDIR/metrics_resampled_data.sh 

# ===============================================
# Distribution of cell type TOTAL number of reads
# ===============================================
while read mark 
do 
    echo "$mark"
    EPDIR=${CHPDIR}/files/${mark}
    cat ${EPDIR}/qc/FINAL_${mark}*.numreads > ${EPDIR}/qc/${mark}_read_totals.tsv
done < <( echo "WCE"; cat ${ANNDIR}/histone_mark_list )

# Make plots: 
R --slave -f $BINDIR/plot_reads_distribution.R

