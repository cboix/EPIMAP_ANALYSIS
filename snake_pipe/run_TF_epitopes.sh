#!/bin/bash
# ---------------------------------------
# Clean and then run TFs or histone marks
# NOTE: Do after WCE is done.
# TODO: Should avoid running finished jobs:
# How should we do this with these array jobs?
# ---------------------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

# Files:
export MARKS=${DBDIR}/Annotation/histone_mark_list
export ALLCHIP=${DBDIR}/Annotation/unordered_epitopes_for_processing
cd $SNDIR

# NOTE: First run WCE, ATAC-seq, and DNase-seq
while read epitope
do
    echo $epitope
    export EPITOPE=$epitope
    ${SNDIR}/clean_empty_tagAligns.sh $epitope
    ${SNDIR}/submit_snakemake_parallel.sh $epitope
# done < <( echo "WCE\nDNase-seq\nATAC-seq" )
done < <( echo "DNase-seq\nATAC-seq" )

while read epitope
do
    echo $epitope
    export EPITOPE=$epitope
    ${SNDIR}/clean_empty_tagAligns.sh $epitope
    ${SNDIR}/submit_snakemake_parallel.sh $epitope
done < $MARKS 

# Run 2ndary epitopes: CTCF and assoc + POLR2A + EP300.
while read epitope
do
    echo $epitope
    export EPITOPE=$epitope
    ${SNDIR}/clean_empty_tagAligns.sh $epitope
    ${SNDIR}/submit_snakemake_parallel.sh $epitope
done < <( echo "CTCF\nEP300\nPOLR2A\nSMC3\nRAD21" )

while read epitope
do
    ISHIST=$( grep "$epitope" $MARKS -c )
    if [[ "$ISHIST" == "0" ]]
    then
        echo $epitope
        ${SNDIR}/clean_empty_tagAligns.sh $epitope
        ${SNDIR}/submit_snakemake_parallel.sh $epitope
    fi
done < $ALLCHIP

# Set up options:
# Command:
while read epitope RERUN; do
    export EPITOPE=$epitope
    SCHEDOPT="-t $RERUN -N ${EPITOPE}_snake -o ${DBDIR}/out/ChIP"
    qsubcmd="qsub $SCHEDOPT ${SNDIR}/run_snakemake.sh general.snake all"
    echo $qsubcmd
    bash -c "$qsubcmd"
done < rerun_list2.tsv

# Eval what needs to be done:
# For rhel7 or rhel6, locate the python libraries:
export LD_LIBRARY_PATH=/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_3.4.2/lib:/broad/software/free/Linux/redhat_7_x86_64/pkgs/python_3.5.1/lib:$LD_LIBRARY_PATH
# rm $DBDIR/out/ChIP/test_log
while read epitope
do
    export EPITOPE=$epitope
    echo $EPITOPE
    cd $SNDIR
    # Count number of cell types:
    if [[  "$EPITOPE" == "DNase-seq" ]] || [[ "$EPITOPE" == "ATAC-seq" ]]; then
        export FILEINFO=${DBDIR}/${EPITOPE}/file_links/all_submitted_released/${EPITOPE}_bam.csv
    else
        export FILEINFO=${DBDIR}/ChIP-seq/file_links/all_submitted_released/${EPITOPE}_bam.csv
    fi
    # Set up options:
    echo "Processing epitope: $EPITOPE - has $( wc -l $FILEINFO | awk '{print $1 - 1}' ) files"
    CTNUM=$( awk -vFS="\t" 'NR>1{print $5}' $FILEINFO | sort -u | wc -l - | awk '{print $1}')
    for task in `seq 1 $CTNUM`; do
        export SGE_TASK_ID=$task
        echo "$task ${EPITOPE}"
        echo "${EPITOPE} $task" >> $DBDIR/out/ChIP/test_log
        snakemake --snakefile general.snake -npk --rerun-incomplete --nolock all &>> $DBDIR/out/ChIP/test_log
    done
done < $MARKS
# NOTE: EVAL ALSO DNase-seq and ATAC-seq
