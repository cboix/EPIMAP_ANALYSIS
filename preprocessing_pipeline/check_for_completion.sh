#!/bin/bash
# Check completion status of epitope or 
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]
then
    export TOPDIR=/broad/compbio/cboix/EPIMAP_ANALYSIS
else
    export TOPDIR=$HOME/EPIMAP_ANALYSIS
fi
export DBDIR=${TOPDIR}/db
export BINDIR=${TOPDIR}/bin
export SNDIR=${BINDIR}/preprocessing_pipeline
export ANNDIR=${DBDIR}/Annotation
cd $SNDIR
TODAY=$(date +%Y%m%d-%H%M)
TMPFILE=${TMP}/${RANDOM}_check_completion_${TODAY}

# Run either all or just one epitope:
if [[ $# -lt 1 ]]
then
    cat ${ANNDIR}/unordered_epitopes_for_processing > $TMPFILE
    CSTAT=${ANNDIR}/completion_status_${TODAY}
    RERUNS=${ANNDIR}/rerun_files_${TODAY}
else
    echo $1 > $TMPFILE
    CSTAT=${ANNDIR}/completion_status_${TODAY}_${1}
    RERUNS=${ANNDIR}/rerun_files_${TODAY}_${1}
fi

while read epitope
do
    echo $epitope
    export EPITOPE=$epitope
    if [[  "$EPITOPE" == "DNase-seq" ]] 
    then
        export FILEINFO=${DBDIR}/${EPITOPE}/file_links/all_submitted_released/${EPITOPE}_bam.tsv
    else
        export FILEINFO=${DBDIR}/ChIP-seq/file_links/all_submitted_released/${EPITOPE}_bam.csv
    fi
    CTNUM=$( awk -vFS="\t" 'NR>1{print $5}' $FILEINFO | sort -u | wc -l - | awk '{print $1}')
    # For each job, check completion:
    for ctid in `seq 1 $CTNUM`
    do 
        export SGE_TASK_ID=$ctid
        if [[ "$EPITOPE" == "WCE" ]] 
        then
            RUNSTAT=$( $SNDIR/run_snakemake.sh --snakefile general.snake -npk --nolock process_pool | grep "Nothing to be done" -c )
        else
            RUNSTAT=$( $SNDIR/run_snakemake.sh --snakefile general.snake -npk --nolock all | grep "Nothing to be done" -c )
        fi
        if (( $RUNSTAT != 1 )) 
        then
            echo -e "$EPITOPE\t$ctid" >> $RERUNS
        fi
        echo -e "$EPITOPE\t$ctid\t$RUNSTAT" >> $CSTAT
    done 
done < $TMPFILE
rm $TMPFILE
popd

# ========================================
# If we want to rerun, FORCING all reruns:
# ========================================
if [[ "0" == "1" ]]
then
    while read epitope id
    do
        export EPITOPE=$epitope
        export FILEINFO=${DBDIR}/ChIP-seq/file_links/all_submitted_released/${EPITOPE}_bam.csv
        # Hold:
        JOBS=$( qstat | wc -l )
        while (( $JOBS > 100 )); do
            sleep 60; JOBS=$( qstat | wc -l )
        done
        # Run:
        SCHEDOPT="-l h_vmem=50G -l h_rt=45:00:00 -P compbio_lab -j y -b y -V -r y -N ${EPITOPE}_snake"
        qsub -cwd -t $id ${SCHEDOPT} -hold_jid WCE_snake -o $DBDIR/out/ChIP "$SNDIR/run_snakemake.sh --snakefile general.snake -pk --nolock all"
    done < $RERUNS
fi

