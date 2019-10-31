#!/bin/bash
if [[ $# -gt 0 ]]; then
    INFOFILE=$1
else
    INFOFILE=${FULLIMPUTATION_TAB}
fi

if [[ $# -gt 1 ]]; then
    SAMPLETABLE=$2
else
    SAMPLETABLE=${SAMPLEMARK_TAB}
fi

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

cd ${DISTANCE_DIR}
wc *.txt -l | awk -vOFS="\t" '$1 == 0{gsub(".txt","",$0); print $2}'  > ${TMP}/tmp_reruns
awk -vOFS="\t" '{print NR,$1"_"$2}' $SAMPLETABLE | join -1 2 -2 1 - ${TMP}/tmp_reruns | awk '{print $2}' > ${TMP}/tmp_rerun_tasks
NRERUN=$( wc -l ${TMP}/tmp_rerun_tasks | awk '{print $1}' )
echo "[STATUS] $NRERUN distance tasks need to be rerun"


if [[ "$NRERUN" != "0" ]];then
    while read item
    do
        qsub -cwd -t $item -P compbio_lab -l h_vmem=8G -l h_rt=8:00:00 -N dist_catch -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh DISTANCE $INFOFILE $SAMPLETABLE"
    done < ${TMP}/tmp_rerun_tasks
fi


# #
# cd ${IMPDIST_DIR}
# wc *.txt -l | awk -vOFS="\t" '$1 == 0{gsub(".txt","",$0); print $2}'  > ${TMP}/tmp_reruns
# awk -vOFS="\t" '{print NR,$1"_"$2}' $IMPOBS_TAB | join -1 2 -2 1 - ${TMP}/tmp_reruns | awk '{print $2}' > ${TMP}/tmp_rerun_tasks
# NRERUN=$( wc -l ${TMP}/tmp_rerun_tasks | awk '{print $1}' )
# echo "[STATUS] $NRERUN distance tasks need to be rerun"


# if [[ "$NRERUN" != "0" ]];then
#     while read item
#     do
#         qsub -cwd -t $item -P compbio_lab -l h_vmem=8G -l h_rt=8:00:00 -N dist_catch -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_obsimp.sh DISTANCE"

# qsub -cwd -t $item -P compbio_lab -tc 2500 -l h_vmem=5G -l h_rt=9:00:00 -N imp_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_obsimp.sh"

#     done < ${TMP}/tmp_rerun_tasks
# fi

