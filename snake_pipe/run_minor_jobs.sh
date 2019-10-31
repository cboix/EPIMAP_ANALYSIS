#!/bin/bash
# ============
# FOR GENERAL: 
# ============
# Run minor jobs
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]
then
    export TOPDIR=/broad/compbio/cboix/EPIMAP_ANALYSIS
else
    export TOPDIR=$HOME/EPIMAP_ANALYSIS
fi
export DBDIR=${TOPDIR}/db
export BINDIR=${TOPDIR}/bin
export SNDIR=${BINDIR}/snake_pipe
export ANNDIR=${DBDIR}/Annotation
TODAY=$(date +%Y%m%d-%H%M)
LOGFILE=${ANNDIR}/minor_jobs_log_${TODAY}

# Run either all or just one epitope:
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
            NJOBS=$( $SNDIR/run_snakemake.sh --snakefile general.snake -npk --nolock process_pool | tail -n 1 | awk '{print $1}' )
        else
            NJOBS=$( $SNDIR/run_snakemake.sh --snakefile general.snake -npk --nolock all | tail -n 1 | awk '{print $1}' )
        fi
        msg="$epitope $ctid has $NJOBS left"
        echo "$msg"
        echo "$msg" >> $LOGFILE
        if (( $NJOBS < 5 ))
        then
            if [[ "$EPITOPE" == "WCE" ]] 
            then
                $SNDIR/run_snakemake.sh --snakefile general.snake -pk --nolock process_pool
            else
                $SNDIR/run_snakemake.sh --snakefile general.snake -pk --nolock all
            fi
        fi
    done 
done < ${TOPDIR}/db/Annotation/unordered_epitopes_for_processing


# Run from list in OUTDIR:
cd $SNDIR
while read epitope id
do
    echo "$id for $epitope"
    export EPITOPE=$epitope
    SCHEDOPT="-N ${EPITOPE}_snake -o ${DBDIR}/out/ChIP"
    if [[  "$EPITOPE" == "DNase-seq" ]] 
    then
        export FILEINFO=${DBDIR}/${EPITOPE}/file_links/all_submitted_released/${EPITOPE}_bam.tsv
    else
        export FILEINFO=${DBDIR}/ChIP-seq/file_links/all_submitted_released/${EPITOPE}_bam.csv
    fi
    # Select run mode:
    if [[ "$EPITOPE" == "WCE" ]] 
    then
        runmode="process_pool"
    elif [[ "$EPITOPE" == "DNase-seq" ]] 
    then
        runmode="all"
    else
        runmode="all"
        SCHEDOPT="$SCHEDOPT -hold_jid WCE_snake"
    fi
    # Command:
    qsubcmd="qsub -t $id $SCHEDOPT ${SNDIR}/run_snakemake.sh general.snake $runmode"
    echo "$qsubcmd"
    bash -c "$qsubcmd"
done < ${TOPDIR}/db/out/ChIP/to_rerun_all

# ==============
# FOR PROJECTS: 
# ==============
export PROJECT=CRs
export EXPT=ChIP-seq
# Directories:
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]
then
    export DBDIR=/broad/compbio/cboix/${PROJECT}/db
else
    export DBDIR=$HOME/${PROJECT}/db
fi
export SFTDIR=/broad/compbio/cboix/software/bin
export BINDIR=$HOME/${PROJECT}/bin
export ANNDIR=$DBDIR/Annotation
export ENCDIR=${HOME}/data/EPIMAP_ANALYSIS
export SNDIR=${ENCDIR}/bin/snake_pipe
export MAPPING=${ENCDIR}/db/Annotation/all_submitted_released_biosample_mapping.tsv
mkdir -p $DBDIR/out $ANNDIR $BINDIR
cd $SNDIR
TODAY=$(date +%Y%m%d-%H%M)
LOGFILE=${ANNDIR}/minor_jobs_log_${TODAY}

export FILEINFO=${ANNDIR}/${PROJECT}_${EXPT}_bam.tsv
echo "Processing project: ${PROJECT} - ${EXPT} has $( wc -l $FILEINFO | awk '{print $1 - 1}' ) files"

CTNUM=$( awk -vFS="\t" -vOFS="\t" 'NR>1{print $3,$4}' $FILEINFO | sort -u | wc -l - | awk '{print $1}')
for ctid in `seq 1 $CTNUM`
do 
    export SGE_TASK_ID=$ctid
    # NJOBS=$( $SNDIR/run_snakemake.sh --snakefile project.snake -npk --nolock all | tail -n 1 | awk '{print $1}' )
    NJOBS=$( snakemake --snakefile project.snake -npk --nolock all | tail -n 1 | awk '{print $1}' )
    msg="$PROJECT $ctid has $NJOBS left"
    echo "$msg"
    echo "$msg" >> $LOGFILE
    if [[ "$NJOBS" != "Nothing" ]]
    then
        if (( $NJOBS < 25 ))
        then
            snakemake --snakefile project.snake -pk --nolock all
            # $SNDIR/run_snakemake.sh --snakefile project.snake -pk --nolock all
        fi
    fi
done 
