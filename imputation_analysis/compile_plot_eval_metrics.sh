#!/bin/bash
# ------------------------------------------
# A. Compile and plot QC metrics of imputation:
# B. Develop metric for quality of imputation 
#       based on observed imputation
# ------------------------------------------
# 0. Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh
TMP_DIR=${TMP}/compile_eval_${RANDOM}
mkdir -p ${TMP_DIR}
cd ${EVAL_DIR}

# 1. Compile all eval:
echo "[STATUS] Aggregating qc metrics"
EVALFILE=${TMP_DIR}/all_eval_${TODAY}.tsv
rm ${EVALFILE}
for file in `ls *_eval.txt`; do
    basename=${file%_*}
    cell=${basename%_*}
    mark=${basename#*_}
    if [[ ! -s $EVALFILE ]]; then
        awk -vOFS="\t" 'NR==1{print "cell","mark",$0}' $file >> ${EVALFILE}
    fi
    awk -vOFS="\t" -v pref="${cell}\t$mark" 'NR==2{print pref,$0}' $file >> ${EVALFILE}
done

# Get the qc metrics for all tracks:
QCFILE=${TMP_DIR}/all_qc_${TODAY}.tsv
rm $QCFILE
while read mark; do
    echo $mark
    if [[ "$mark" == "DNase-seq" ]] || [[ "$mark" == "ATAC-seq" ]]; then
        QCDIR=$DBDIR/${mark}/files/${mark}/qc
    else
        QCDIR=$DBDIR/ChIP-seq/files/${mark}/qc
    fi
    cat $QCDIR/FINAL_${mark}_*sub.cc.qc >> ${QCFILE}
done < ${MARKS_LIST}

# 2. Plot eval metrics:
# (1) Table of eval metric rank: TODO
# (2) Barplots of metrics
echo "[STATUS] Plotting qc metrics"
source activate mv_env
R --slave -f ${BINDIR}/plot_chromImpute_metrics.R --args ${EVALFILE} ${QCFILE}
source deactivate

# (3) Plot obs/imputed tracks (top/bottom 5) for each mark:
source activate mv_env > /dev/null
for mark in `awk 'NR > 1{print $2}' $EVALFILE | sort -u`; do
    echo $mark
    awk -v mark=$mark 'NR > 1 && $2 == mark {print $1, $2, $6}' $EVALFILE | sort -k3 -n > ${TMP_DIR}/${mark}_tmp
    MARKFILE=${TMP_DIR}/${mark}_top_bottom_cor.tsv
    head ${TMP_DIR}/${mark}_tmp -n 5 | awk '{print $1"_"$2, NR}' > $MARKFILE
    tail ${TMP_DIR}/${mark}_tmp -n 5 | awk '{print $1"_"$2, NR + 5}' >> $MARKFILE
    sort -u $MARKFILE > ${MARKFILE}.tmp
    awk -vOFS="\t" '{print $1"_"$2, $3}' $SAMPLEMARK_TAB | sort -u | join - ${MARKFILE}.tmp | sort -k3 -n | awk -vOFS="\t" '{gsub("_","\t", $1); print $1,$2}' > $MARKFILE

    RPATH="/broad/compbio/cboix/software/miniconda2/envs/mv_env/bin/R"
    qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=1:30:00 -N impobs_$mark -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; $RPATH --slave -f ${BINDIR}/plot_impobs_tracks.R --args ${MARKFILE} chr10 0 135000000"
    qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=1:30:00 -N impobs_$mark -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; $RPATH --slave -f ${BINDIR}/plot_impobs_tracks.R --args ${MARKFILE} chr1 1500000 2000000"
done
source deactivate > /dev/null


# TODO: Other metrics sources (like enrichment) should be plotted.

end=`date +%s`
runtime=$((end-start))
echo "Finished enrichment pipeline in $runtime seconds."
