#!/bin/bash
# ====================================================
# Process DNase and ChIP-seq epitopes using snakemake:
# (NOTE: Should use process large batch epitopes only)
# ====================================================
export EPITOPE=$1
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 15

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
SCHEDOPT="-t 1-$CTNUM -N ${EPITOPE}_snake -o ${DBDIR}/out/ChIP"
# Choose run mode:
if [[ "$EPITOPE" == "WCE" ]]; then
    runmode="process_pool"
elif [[  "$EPITOPE" == "DNase-seq" ]] || [[ "$EPITOPE" == "ATAC-seq" ]]; then
    runmode="all"
else
    runmode="all"
    SCHEDOPT="$SCHEDOPT -hold_jid WCE_snake"
fi

# Command:
qsubcmd="qsub $SCHEDOPT ${SNDIR}/run_snakemake.sh general.snake $runmode"

echo "[STATUS] Running the qsub array job:" 
echo "$qsubcmd"
bash -c "$qsubcmd"

# # Eval what needs to be done:
for task in `seq 1 $CTNUM`; do
    export SGE_TASK_ID=$task
    echo $task
    snakemake --snakefile general.snake -pk --rerun-incomplete --nolock all
done

# # Eval what needs to be done:
# for task in `seq 221 $CTNUM`; do
#     export SGE_TASK_ID=$task
#     echo $task
#     snakemake --snakefile general.snake -pk --nolock process_pool 
#     #snakemake --snakefile general.snake -npk --nolock process_pool | tail -n 4 >> ${ANNDIR}/taillogs_${EPITOPE}_${TODAY}
# done
