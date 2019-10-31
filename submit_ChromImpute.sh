#!/bin/bash
# ======================================================
# Run ChromImpute in parallel chunks from converted data
# ChromImpute code by J.Ernst:
# http://www.biolchem.ucla.edu/labs/ernst/ChromImpute/
# ======================================================
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

RUNACC=0  # RUN CTCF/ETC.
RUNACC=1  # RUN CTCF/ETC.

# ===================================================
# Make table of existing files in converted data dir:
# NOTE: Fine to recreate this, as it takes stock of 
#       whatever is currently available.
# NOTE: Conversion is done in the snakemake pipeline.
# ===================================================
# Format as: sample1 mark1 fileA
if [[ ! -s $SAMPLEMARK_TAB ]]; then
    ls ${CONVERTED_DIR}/chr1_FINAL*BSS[0-9]*pval.signal.bedgraph.gz.wig.gz | awk -v OFS="\t" -v cvdir="${CONVERTED_DIR}/" '{a = $0; sub(cvdir, "", a); sub("chr[0-9XY]*_FINAL_","",a);
    sub("\.wig\.gz$","",a);
    cell_mark = a; sub ("\.sub_VS.*","",cell_mark);
    cell = cell_mark; sub("[0-9A-Za-z-]*_","",cell);
    mark = cell_mark; sub("_.*$","",mark); print cell, mark ,"FINAL_"a}' | sort -u > ${SAMPLEMARK_TAB}

    # EXTENDED with CTCF, etc:
    awk '$2 ~ /^H[1234]/||/CTCF/||/POLR2A/||/SMC3/||/RAD21/||/EP300/||/DNase-seq/||/ATAC-seq/' ${SAMPLEMARK_TAB} > ${SAMPLEMARK_TAB}_ext  # What about DNase...

    # Cut down file to just histone modifications:
    # NOTE: keeps CTCFL - we can see what that looks like?
    awk '$2 ~ /^H[1234]/||/DNase-seq/||/ATAC-seq/' ${SAMPLEMARK_TAB} > ${SAMPLEMARK_TAB}_tmp
    mv ${SAMPLEMARK_TAB}_tmp ${SAMPLEMARK_TAB}
    awk '{print $1}' ${SAMPLEMARK_TAB} | sort -u > ${CELLINFO}

    # Prune ext file:
    sort +0 -1 ${SAMPLEMARK_TAB}_ext | join ${CELLINFO} - | awk -vOFS="\t" '{print $1,$2,$3}' > ${EXT_SAMPMARK_TAB}
    rm ${SAMPLEMARK_TAB}_ext
fi

# Impute 13 top marks OR impute necessary marks for chromHMM
if [[ ! -s ${TOIMPUTE_TAB} ]] || [[ ! -s ${FULLIMPUTATION_TAB} ]];then
    rm -rf ${TOIMPUTE_TAB} ${FULLIMPUTATION_TAB}
    cut -f2 ${SAMPLEMARK_TAB} | sort | uniq -c | sort | tail -n 13 | awk '{print $2}' > ${CIDIR}/top_marks.tsv
    while read celltype
    do
        for mark in ${epitopes[@]}; do
            line="$celltype\t$mark"
            if ! grep -P "$line" ${SAMPLEMARK_TAB} > /dev/null; then
                echo "$line" >> ${TOIMPUTE_TAB}
            fi
            echo "$line" >> ${FULLIMPUTATION_TAB}
        done
    done < <( cut -f1 ${SAMPLEMARK_TAB} | sort -u )

    # Extra 7 marks:
    rm -rf ${EXTRAIMPUTATION_TAB} 
    cut -f2 ${SAMPLEMARK_TAB} | sort | uniq -c | sort -r | head -n 13 | awk '{print $2}' > ${CIDIR}/top_marks.tsv
    while read mark; do
        if [[ "$( grep $mark $FULLIMPUTATION_TAB | wc -l | awk '{print $1}' )" == "0"  ]]; then
            while read celltype
            do
                line="$celltype\t$mark"
                echo "$line" >> ${EXTRAIMPUTATION_TAB}
            done < <( cut -f1 ${SAMPLEMARK_TAB} | sort -u )
        fi
    done < ${CIDIR}/top_marks.tsv
    sort $EXTRAIMPUTATION_TAB > $EXTRAIMPUTATION_TAB.tmp
    mv $EXTRAIMPUTATION_TAB.tmp $EXTRAIMPUTATION_TAB 

    # Accessory:
    rm -rf ${ACCSIMPUTATION_TAB} 
    for mark in ${accs_epitopes[@]}; do
        while read celltype
        do
            line="$celltype\t$mark"
            echo "$line" >> ${ACCSIMPUTATION_TAB}
        done < <( cut -f1 ${SAMPLEMARK_TAB} | sort -u )
    done
    sort $ACCSIMPUTATION_TAB > ${ACCSIMPUTATION_TAB}.tmp
    mv ${ACCSIMPUTATION_TAB}.tmp $ACCSIMPUTATION_TAB
fi
INUM=$( wc -l ${TOIMPUTE_TAB} | awk '{print $1}' )

# NOTE: Run separately
if [[ "$NUMSTATES" == "25" ]]; then 
    export INFOFILE=$EXTRAIMPUTATION_TAB
    export SAMPLETABLE=${SAMPLEMARK_TAB}
else
    if [[ "$RUNACC" == "1" ]]; then
        export INFOFILE=$ACCSIMPUTATION_TAB
        export SAMPLETABLE=${EXT_SAMPMARK_TAB}
    else
        export INFOFILE=$FULLIMPUTATION_TAB
        export SAMPLETABLE=${SAMPLEMARK_TAB}
    fi
fi
FNUM=$( wc -l ${INFOFILE} | awk '{print $1}' )
echo "New (Core): $INUM, torun: $FNUM"

CTNUM=$( cut -f1 $SAMPLETABLE | sort -u | wc -l)
SNUM=$( wc -l ${SAMPLETABLE} | awk '{print $1}' )
awk -v FS="\t" '{print $2}' ${SAMPLETABLE} | sort -u > ${MARKS_LIST}
NMARKS=$( wc -l ${MARKS_LIST} | awk '{print $1}' )
echo "$SNUM files for $NMARKS marks in $CTNUM biosamples"

# ==============================================
# 2. Compute global distance between datasets
# NOTE: For each sample + mark, ranks nearest samples by correlation of the mark.
#       Can parallelize by both sample and mark or by just mark: 
#       flags: -s $SAMPLE $MARK or just -m $MARK
# ==============================================
# Takes ~ 3:30 on slower computers for DNase-seq.
qsub -cwd -t 1-$SNUM -P compbio_lab -tc 2500 -l h_vmem=5G -l h_rt=7:00:00 -N dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh DISTANCE ${INFOFILE} ${SAMPLETABLE}"

qsub -cwd -P compbio_lab -l h_vmem=5G -l h_rt=0:30:00 -N dist_catch -j y -b y -V -r y -o $DBDIR/out/ChromImpute -hold_jid dist_${UQID} "${BINDIR}/catch_rerun.sh $INFOFILE $SAMPLEMARK_TAB"

# FOR RERUNS:
if [[ "1" == "0" ]];then
    cd ${DISTANCE_DIR}
    wc *.txt -l | awk -vOFS="\t" '$1 == 0{gsub(".txt","",$0); print $2}'  > ${TMP}/tmp_reruns
    awk -vOFS="\t" '{print NR,$1"_"$2}' $SAMPLEMARK_TAB | join -1 2 -2 1 - ${TMP}/tmp_reruns | awk '{print $2}' > ${TMP}/tmp_rerun_tasks
    NRERUN=$( wc -l ${TMP}/tmp_rerun_tasks | awk '{print $1}' )
    echo "[STATUS] $NRERUN distance tasks need to be rerun"


    if [[ "$NRERUN" != "0" ]];then
        while read item
        do
            qsub -cwd -t $item -P compbio_lab -l h_vmem=8G -l h_rt=8:00:00 -N dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh DISTANCE"
        done < ${TMP}/tmp_rerun_tasks
    fi
fi

# =====================================================
# 3. Generate training features (separately for each mark)
# NOTE: Given data and correlations, generates training data instances.
#       Can parallelize by using -c $CHROM over each chromsome.
#       Is already broken up by marks.
#       Can specify # loc for training: -f (default 100,000)
#       For training: set fixed -d $SEED 
# =====================================================
# export CHROM=""
# qsub -cwd -t 1-$NMARKS -P compbio_lab -tc 2500 -l h_vmem=40G -l h_rt=15:00:00 -N feat_${UQID} -j y -b y -V -r y -hold_jid dist_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh FEATURES"

# SHORT BY CHROMOSOME:
IFS=$'\t'
while read chr size
do
    export CHROM=$chr
    echo $chr
    # NOTE: Name is not unique so that pred_${UQID} has to wait for all to finish.
    # 5 hrs is fine, except for a few runs that take a long time and truncate...
    qsub -cwd -t 1-$NMARKS -P compbio_lab -tc 2500 -l h_vmem=30G -l h_rt=8:00:00 -N feat_${UQID} -j y -b y -V -r y -hold_jid dist_catch -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh FEATURES $INFOFILE $SAMPLETABLE"
done < ${CHROMSIZES}
export CHROM=""
 
# NOTE: Make sure to rerun all that do not end. Check for EOF errors in output files.

# ===============================================================
# 4. Generate the trained predictors for a specific mark + sample
# NOTE: Trains regression trees based on training data
#       Parallelization enforced by sample/mark combinations.
#       Must match options to GenerateTrainData
# ===============================================================
# NOTE: 4 hours is enough - Takes about 20 min-1hr on ish, but some cores are terribly slow...
qsub -cwd -t 1-$FNUM -P compbio_lab -tc 2500 -l h_vmem=8G -l h_rt=07:00:00 -N pred_${UQID} -j y -b y -V -r y -hold_jid feat_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh PREDICTORS $INFOFILE $SAMPLETABLE"

# For reruns:
if [[ "1" == "0" ]]; then
    while read item
    do
        qsub -cwd -t $item -P compbio_lab -l h_vmem=10G -l h_rt=08:00:00 -N pred_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh PREDICTORS $INFOFILE $SAMPLETABLE"
    done < failed_pred_runs
fi

# ========================================================
# 5. Generate the imputed signal track for a mark + sample
# NOTE: Parallelized over sample/marks
#       Can parallelize over chromosomes with -c $CHROM
# ========================================================
# Impute full set of chrom::
# For all chrom, needed ~48 hrs
# qsub -cwd -t 1-$FNUM -P compbio_lab -l h_vmem=10G -l h_rt=48:00:00 -N impute_${UQID} -j y -b y -V -r y -hold_jid_ad pred_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh IMPUTE"
IFS=$'\t'
while read chr size
do
    export CHROM=$chr
    # 5 hrs should be enough - unless chr1 takes very long. 
    # Matched to predictors array job + can be in parallel
    qsub -cwd -t 1-$FNUM -P compbio_lab -l h_vmem=9G -l h_rt=7:30:00 -N impute_${UQID} -j y -b y -V -r y -hold_jid_ad pred_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh IMPUTE $INFOFILE $SAMPLETABLE"
done < ${CHROMSIZES}
export CHROM="" # Must clear chromosome value.
# NOTE: Make sure to rerun all that do not end. Check for EOF errors in output files.

# =======================================================
# 6. Evaluate imputation performance for existing tracks:
# =======================================================
qsub -cwd -t 1-$SNUM -P compbio_lab -tc 2500 -l h_vmem=5G -l h_rt=2:45:00 -N eval_${UQID} -hold_jid impute_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh EVAL $INFOFILE $SAMPLETABLE"

# =================================================================
# 7. Turn Imputed tracks into binarized for ChromHMM:
# =================================================================
qsub -cwd -t 1-$FNUM -P compbio_lab -tc 2500 -l h_vmem=10G -l h_rt=00:30:00 -N export_${UQID} -j y -b y -V -r y -hold_jid_ad impute_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh EXPORT $INFOFILE $SAMPLETABLE"

# Export also with new cutoffs:
qsub -cwd -t 1-$FNUM -P compbio_lab -tc 2500 -l h_vmem=8G -l h_rt=00:40:00 -N qcut_export_${UQID} -j y -b y -V -r y -hold_jid_ad impute_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh QEXPORT $INFOFILE $SAMPLETABLE"

# 7a. Look at distribution for the tracks 
UNUM=$( wc -l ${ALL_UQ_TAB} | awk '{print $1}' )
qsub -cwd -t 1-$UNUM -P compbio_lab -l h_vmem=3G -l h_rt=00:30:00 -N getsignaldistr_${UQID} -j y -b y -V -r y -hold_jid_ad impute_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/run_eval_signal_distribution.sh $ALL_UQ_TAB $DISTR_DIR"


# ===========================================
# 8. Run ChromHMM for mixed expt and imputed:
# ===========================================
# NOTE: Make impobs and mixobs tables first.
source $BINDIR/make_impobs_tables.sh
MODELS=(15 18 25)
for NUMSTATES in ${MODELS[@]}; do
    echo $NUMSTATES
    # Only imputed
    # CELLNUM=$( cut -f1 $IMPOBS_TAB | sort -u | wc | awk '{print $1}' )
    # qsub -cwd -t 1-$CELLNUM -P compbio_lab -tc 2500 -l h_vmem=10G -l h_rt=00:45:00 -N annotate_${NUMSTATES}_${UQID} -j y -b y -V -r y -hold_jid export_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromHMM.sh -i ${IMPOBS_TAB} -n $NUMSTATES"
    # Only imputed, with binarization according to avg signal distr.:
    # CELLNUM=$( cut -f1 $IMPOBS_TAB | sort -u | wc | awk '{print $1}' )
    # qsub -cwd -t 1-$CELLNUM -P compbio_lab -tc 2500 -l h_vmem=10G -l h_rt=00:45:00 -N annotate_qcut_${NUMSTATES}_${UQID} -j y -b y -V -r y -hold_jid qcut_export_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromHMM.sh -i ${IMPOBS_TAB} -n $NUMSTATES -q 1"

    # Mixed imputed and observed:
    MIXNUM=$( cut -f1 $MIXOBS_TAB | sort -u | wc | awk '{print $1}' )
    # qsub -cwd -t 1-$CELLNUM -P compbio_lab -tc 2500 -l h_vmem=10G -l h_rt=00:45:00 -N annotate_${NUMSTATES}_${UQID} -j y -b y -V -r y -hold_jid export_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromHMM.sh -i ${MIXOBS_TAB} -n $NUMSTATES -p mixed_impobs"
    qsub -cwd -t 1-$MIXNUM -P compbio_lab -tc 2500 -l h_vmem=10G -l h_rt=00:45:00 -N annotate_qcut_mix_${NUMSTATES}_${UQID} -j y -b y -V -r y -hold_jid export_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromHMM.sh -i ${MIXOBS_TAB} -n $NUMSTATES -p mixed_impobs -q 1"
done

# ===========================================
# 9b. Run ChromHMM with new models
# ===========================================

# ==================================
# 9. Call peaks from imputed tracks:
# ==================================
qsub -cwd -t 1-$FNUM -P compbio_lab -tc 2500 -l h_vmem=10G -l h_rt=01:45:00 -N callpkimp_${UQID} -j y -b y -V -r y -hold_jid_ad impute_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/callpeak_from_imputedbw.sh $INFOFILE $PKIMP_DIR"

EPITOPE=DNase-seq
awk -v epitope=$EPITOPE '$2 == epitope' $INFOFILE > $INFOFILE.tmp
DNUM=$( awk 'END{print NR}' $INFOFILE.tmp )
echo "[STATUS] Running call pk impute for $DNUM $EPITOPE tracks"
qsub -cwd -t 1-$DNUM -P compbio_lab -tc 2500 -l h_vmem=10G -l h_rt=01:45:00 -N callpkimp_${EPITOPE}_${UQID} -j y -b y -V -r y -hold_jid_ad impute_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/callpeak_from_imputedbw.sh $INFOFILE.tmp $PKIMP_DIR"

# ===================================
# 10. Call peaks from observed tracks:
# ===================================
qsub -cwd -t 1-$SNUM -P compbio_lab -tc 2500 -l h_vmem=10G -l h_rt=01:30:00 -N callpkobs_${UQID} -j y -b y -V -r y -hold_jid_ad impute_${UQID} -o $DBDIR/out/ChromImpute "$BINDIR/callpeak_from_imputedbw.sh $SAMPLEMARK_TAB $PKOBS_DIR ${CONVERTED_DIR}"

# ==================================
# 11. Plot some areas of the genome:
# ==================================
MODNAMES=(imputed12_25 imputed12_25_on_mixed_impobs observed_15 observed_15_on_mixed_impobs observed_aux_18 observed_aux_18_on_mixed_impobs)
for MOD in ${MODNAMES[@]}; do
    echo $MOD
    base=${MOD%%_on_*}
    NSTATES=${base##*_}
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr1 -a 1000000 -b 1250000 -n $NSTATES -m $MOD"
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr1 -a 1000000 -b 2000000 -n $NSTATES -m $MOD"
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr2 -a 0 -b 5000000 -n $NSTATES -m $MOD"
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr10 -a 0 -b 2000000 -n $NSTATES -m $MOD"
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr19 -a 2000000 -b 4000000 -n $NSTATES -m $MOD"
done

# Plot the models:
MODNAMES=(imputed12_25_QCUT observed_15_QCUT observed_aux_18_QCUT imputed12_25 observed_15 observed_aux_18 imputed12_25_on_mixed_impobs observed_15_on_mixed_impobs observed_aux_18_on_mixed_impobs)
for MOD in ${MODNAMES[@]}; do
    base=${MOD%%_QCUT}
    base=${base%%_on_*}
    NSTATES=${base##*_}
    echo "$MOD in $NSTATES states"
    # Different loc / resolutions:
    qsub -cwd -P compbio_lab -l h_vmem=8G -l h_rt=1:00:00 -N plotchmm_${MOD} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr19 -a 200000 -b 500000 -n $NSTATES -m $MOD"
    qsub -cwd -P compbio_lab -l h_vmem=8G -l h_rt=1:00:00 -N plotchmm_${MOD} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr1 -a 1000000 -b 1250000 -n $NSTATES -m $MOD"
    qsub -cwd -P compbio_lab -l h_vmem=8G -l h_rt=1:00:00 -N plotchmm_${MOD} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr19 -a 1000000 -b 1500000 -n $NSTATES -m $MOD"
    qsub -cwd -P compbio_lab -l h_vmem=8G -l h_rt=1:00:00 -N plotchmm_${MOD} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr10 -a 1000000 -b 3000000 -n $NSTATES -m $MOD"
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr1 -a 1000000 -b 2000000 -n $NSTATES -m $MOD"
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr2 -a 0 -b 5000000 -n $NSTATES -m $MOD"
done

# LiftOver models:
for MOD in ${MODNAMES[@]}; do
    # NOTE: CHANGE TABLE BASED ON TYPE OF RUN.
    MOD=observed_aux_18_on_mixed_impobs_QCUT
    base=${MOD%%_QCUT}
    base=${base%%_on_*}
    NSTATES=${base##*_}
    echo "$MOD in $NSTATES states"

    MIXNUM=$( cut -f1 $MIXOBS_TAB | sort -u | wc | awk '{print $1}' )
    qsub -cwd -t 1-$MIXNUM -P compbio_lab -l h_vmem=10G -l h_rt=0:30:00 -N lo_chmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_liftOver_on_ChromHMM.sh -n $NSTATES -m $MOD -i $MIXOBS_TAB"
done

# Plot epilogos salient regions (with Roadmap):
while read chr start end; do
    MOD=observed_aux_18_on_mixed_impobs_QCUT
    base=${MOD%%_QCUT}
    base=${base%%_on_*}
    NSTATES=${base##*_}
    echo "$MOD in $NSTATES states $chr $start $end"
    # Different loc / resolutions:
    qsub -cwd -P compbio_lab -l h_vmem=8G -l h_rt=1:00:00 -N plotchmm_${MOD} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c $chr -a $start -b $end -n $NSTATES -m $MOD"
done < $DBDIR/epilogos_regions_101019.tsv

# =====================================================
# 12. Merge masterlist for DNase-seq (and other marks):
# =====================================================
EPITOPE="DNase-seq"
# EPITOPE="H3K27ac"
# EPITOPE="H3K4me1"
# EPITOPE="H3K4me3"

qsub -cwd -P compbio_lab -l h_vmem=30G -l h_rt=3:00:00 -N run_ML_${EPITOPE} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_Index_masterlist_gen.sh -m ${EPITOPE}"

# -------------------------------------------------
# 13. Compress BW files into HDF5 files.
# NOTE: re-run faster the 2nd time: CP concurrently
# -------------------------------------------------
IFS=$'\t'
while read epitope; do 
    while read chr size;do
        # Imputed 
        qsub -cwd -P compbio_lab -l h_vmem=30G -l h_rt=2:00:00 -N comp_imp_${epitope}_${chr} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate pytorch_env > /dev/null; python $BINDIR/compress_bw_to_hd5py.py main --chrom ${chr} --mark ${epitope} --dataset imputed"

        # Observed
        qsub -cwd -P compbio_lab -l h_vmem=30G -l h_rt=2:00:00 -N comp_obs_${epitope}_${chr} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate pytorch_env > /dev/null; python $BINDIR/compress_bw_to_hd5py.py main --chrom ${chr} --mark ${epitope} --dataset observed"
    done < ${CHROMSIZES_noY}
done < <( echo "DNase-seq\nH3K27ac\nH3K27me3\nH3K36me3\nH3K4me1\nH3K4me3\nH3K9me3" )

