#!/bin/bash
# ======================================================
# Run ChromImpute in parallel chunks from converted data
# ChromImpute code by J.Ernst:
# http://www.biolchem.ucla.edu/labs/ernst/ChromImpute/
# ======================================================
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

# ===============================================
# A. Make samplesheet of all imputed and observed
# =============================================
source $BINDIR/make_impobs_tables.sh

rm ${ALL_UQ_TAB}.tmp
extra_epitopes=(H3K9ac H3K79me2 H4K20me1 H3K4me2 H2AFZ ATAC-seq)
# Either run main:
# while read epitope; do
#     echo "Adding: $epitope"
#     grep ${epitope} ${ALL_UQ_TAB} >> ${ALL_UQ_TAB}.tmp
# done < <( echo "DNase-seq\nH3K27ac\nH3K27me3\nH3K36me3\nH3K4me1\nH3K4me3\nH3K9me3" )

# Or run extra:
for epitope in ${extra_epitopes[@]}; do
    echo "Adding: $epitope"
    grep ${epitope} ${ALL_UQ_TAB} >> ${ALL_UQ_TAB}.tmp
done

CTNUM=$( cut -f1 ${ALL_UQ_TAB}.tmp | sort -u | wc -l)
SNUM=$( wc -l ${ALL_UQ_TAB}.tmp | awk '{print $1}' )
echo "$SNUM imputed/observed files in $CTNUM biosamples (split imp/obs)"

# ----------------------------------------------
# 3. Compute global distance between datasets
# ----------------------------------------------
# Takes longer than before because there are more cells to compare.
# NOTE: ALL CONVERTED FILES MUST BE COPIED TO IMPUTE DIRECTORY.

qsub -cwd -t 1-$SNUM -P compbio_lab -tc 2500 -l h_vmem=8G -l h_rt=8:00:00 -N imp_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_obsimp.sh ${ALL_UQ_TAB}.tmp"

# -------------------------
# 4. Alternative distances:
# -------------------------
# a. Preprocess data:
UNUM=$( wc -l ${ALL_UQ_TAB} | awk '{print $1}' )
qsub -cwd -t 1-$UNUM -P compbio_lab -l h_vmem=2G -l h_rt=0:20:00 -N preprocess_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_obsimp.sh PREPROCESS ${ALL_UQ_TAB} ${REGDIST_DIR} 20000"

# b. Run distance fixed all bins
qsub -cwd -t 1-$UNUM -P compbio_lab -l h_vmem=4G -l h_rt=0:45:00 -N reg_distall_${UQID} -hold_jid preprocess_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_obsimp.sh FIXEDDIST ${ALL_UQ_TAB} ${REGDIST_DIR} 0"

# c. Run distance fixed 20000 bins
qsub -cwd -t 1-$UNUM -P compbio_lab -l h_vmem=3G -l h_rt=0:30:00 -N reg_dist20k_${UQID} -hold_jid preprocess_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_obsimp.sh FIXEDDIST ${ALL_UQ_TAB} ${REGDIST_DIR} 20000"

# d. Run distance fixed 100000 bins
qsub -cwd -t 1-$UNUM -P compbio_lab -l h_vmem=3G -l h_rt=0:30:00 -N reg_dist100k_${UQID} -hold_jid preprocess_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_obsimp.sh FIXEDDIST ${ALL_UQ_TAB} ${REGDIST_DIR} 100000"

# e. Run distance fixed 250000 bins
qsub -cwd -t 1-$UNUM -P compbio_lab -l h_vmem=3G -l h_rt=0:30:00 -N reg_dist250k_${UQID} -hold_jid preprocess_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_obsimp.sh FIXEDDIST ${ALL_UQ_TAB} ${REGDIST_DIR} 250000"

# g. Run distance fixed 1M bins
qsub -cwd -t 1-$UNUM -P compbio_lab -l h_vmem=3G -l h_rt=0:30:00 -N reg_dist1M_${UQID} -hold_jid preprocess_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_obsimp.sh FIXEDDIST ${ALL_UQ_TAB} ${REGDIST_DIR} 1000000"

# Collate all into dataframes per mark and tarball:
NREGIONS=0
NREGIONS=20000
NREGIONS=100000
NREGIONS=250000
NREGIONS=1000000
cd ${REGDIST_DIR}
while read mark; do
    echo $mark
    ALLTAB=${REGDIST_DIR}/all_fixeddist_${NREGIONS}_${mark}.tsv
    rm $ALLTAB
    while read file; do
        cat $file >> $ALLTAB
    done < <(ls fixeddist_${NREGIONS}_*${mark}.tsv)
    gzip -f $ALLTAB
done < ${MARKS_LIST}
tar -cvzf allmark_fdist_${NREGIONS}.tar.gz all_fixeddist_${NREGIONS}_*.tsv.gz
popd

NREGIONLIST=(0 250000 1000000)
MEASURES=(spearman cor jaccard)
for nregion in ${NREGIONLIST[@]}; do
    for measure in ${MEASURES[@]}; do 
        R --slave -f $BINDIR/plot_chromImpute_imputed_umap.R --args TRUE $nregion $measure
    done
done

# ---------------------------------------------------------
# 4. Plot all of the distances as heatmaps, trees, and MDS:
# ---------------------------------------------------------
# NOTE: See plot_chromImpute_imputed_distances.R
# ------------------------------------
# 5. Compare to original distances,
# - look at most surprising distances.
# ------------------------------------

# ==============================================
# B. Make samplesheet that prioritizes observed 
# at good QC, but fills gaps with imputed.
# ==============================================
AGTAB=${DBDIR}/best_worst_impobs_agreement_011119.tsv
if [[ ! -s $MIXOBS_TAB ]]; then
    # Start with only imputed data:
    ls ${IMPUTED_DIR}/chr1_impute_BSS*.gz | awk -v OFS="\t" -v cvdir="${IMPUTED_DIR}/" '{a = $0; sub(cvdir, "", a); sub("chr[0-9XY]*_impute_","",a); 
    sub("\.wig\.gz$","",a);
    print a,"impute_"a}' | sort -u > ${MIXOBS_TAB}

    # Remove tracks with good obs/imp aggreement from table (QC choice to use observed):
    awk -vOFS="\t" '$3 == "MEDIUM" || $3 == "TOP" {print $1"_"$2, $3}' $AGTAB | sort -u > $AGTAB.tmp 
    join -v 1 $MIXOBS_TAB $AGTAB.tmp | awk -vOFS="\t" '{sub("_","\t", $1);print $1,$2,$3}' > ${MIXOBS_TAB}.tmp

    # Fill in again with observed data:
    awk -vOFS="\t" '{print $1"_"$2, $3}' $SAMPLEMARK_TAB | sort -u > ${SAMPLEMARK_TAB}.tmp
    join $MIXOBS_TAB $AGTAB.tmp | sort -u | join ${SAMPLEMARK_TAB}.tmp - | awk -vOFS="\t" '{sub("_","\t", $1);print $1,$2}' | tee ${CIDIR}/obs_toadd.tsv >> ${MIXOBS_TAB}.tmp

    # Keep only top epitopes:
    rm ${MIXOBS_TAB}
    for epitope in ${epitopes[@]}; do
        echo $epitope
        grep "_${epitope}" ${MIXOBS_TAB}.tmp >> ${MIXOBS_TAB}
    done 

    # 2. Copy these specific unimputed tracks to impute dir:
    cd $CONVERTED_DIR
    echo "[STATUS] Copying over observed files to cover missing impute"
    while read cell mark suffix; do
        echo $cell $mark
        cp $( ls chr*_${suffix}* | grep 'chr[0-9XY]*_FINAL_*' ) ${IMPUTED_DIR}
    done < ${CIDIR}/obs_toadd.tsv

    sort -k1 -k2 $MIXOBS_TAB > ${MIXOBS_TAB}.tmp
    mv ${MIXOBS_TAB}.tmp ${MIXOBS_TAB}
    rm ${SAMPLEMARK_TAB}.tmp
fi

CTNUM=$( cut -f1 $MIXOBS_TAB | sort -u | wc -l)
SNUM=$( wc -l ${MIXOBS_TAB} | awk '{print $1}' )
echo "$SNUM imputed/observed files in $CTNUM biosamples"

# ----------------------------------------------
# 3. Compute global distance between datasets
# ----------------------------------------------
# Takes longer than before because there are more cells to compare.
qsub -cwd -t 1-$SNUM -P compbio_lab -tc 2500 -l h_vmem=5G -l h_rt=9:00:00 -N imp_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_obsimp.sh $MIXOBS_TAB $MIXDIST_DIR"


