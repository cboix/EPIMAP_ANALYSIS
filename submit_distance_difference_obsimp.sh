#!/bin/bash
# --------------------------------------------
# Run distance metrics within an imputed cell 
# for real vs. difference tracks.
# --------------------------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

TMP_DIR=${TMP}/impute_obs_diff_${RANDOM}
mkdir -p ${TMP_DIR} ${DIFFDIST_DIR}
SNUM=$( wc -l ${SAMPLEMARK_TAB} | awk '{print $1}' )
echo "There are $SNUM observed files"

# Merge impobs and sampletab:
awk '$3 ~ /impute/{print $1"_"$2, $3}' ${ALL_TRACKS_TAB} | sort -u > ${TMP_DIR}/tmp_impobstable
awk '{print $1"_"$2, $3}' ${SAMPLEMARK_TAB} | sort -u | join - ${TMP_DIR}/tmp_impobstable | awk -vOFS="\t" '{sub("_","\t", $1); print $1,$2}' > ${TMP_DIR}/imputed_tracks.tsv
# Match to original prefixes:
INUM=$( wc -l ${TMP_DIR}/imputed_tracks.tsv | awk '{print $1}' )
echo "There are $INUM imputed files that match files"

# ---------------------------------------
# 1. Create all observed - imputed tracks
# ---------------------------------------
qsub -cwd -t 1-$INUM -P compbio_lab -tc 2500 -l h_vmem=5G -l h_rt=2:00:00 -N calcdiff_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_difference_obsimp.sh SCALEDIFF ${TMP_DIR}/imputed_tracks.tsv"

# -----------------------------------------------------
# 2. Run the global distance measures within cell types
# -----------------------------------------------------
# Within-mark:
qsub -cwd -t 1-$INUM -P compbio_lab -tc 2500 -l h_vmem=8G -l h_rt=8:00:00 -N diff_dist_${UQID} -hold_jid_ad calcdiff_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_difference_obsimp.sh SAMPDIST ${TMP_DIR}/imputed_tracks.tsv"

# Within-sample
qsub -cwd -t 1-$INUM -P compbio_lab -tc 2500 -l h_vmem=8G -l h_rt=2:30:00 -N diff_dist_${UQID} -hold_jid_ad calcdiff_dist_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_on_difference_obsimp.sh DISTANCE ${TMP_DIR}/imputed_tracks.tsv"


# --------------------------------------------------------
# 3. Find tracks with a residual that strongly corresponds
# to another mark signal (potential bad antibodies)
# NOTE: Only tracks that actually have imputed:
# --------------------------------------------------------
# a. Aggregate table:
echo "[STATUS] Aggregating distances for difference "
DISTFILE=${TMP_DIR}/all_diffdist_${TODAY}.tsv
SDISTFILE=${TMP_DIR}/all_sampdiffdist_${TODAY}.tsv
echo -e "sample\\tmark\\tagainst\\tdist" > ${DISTFILE}
while read sample mark obsfile; do 
    echo "${sample} ${mark}"
    pref=diff_${mark}_${sample}
    DMAT=${DIFFDIST_DIR}/${pref}.txt
    if [[ -s $DMAT ]]; then 
        awk -v sample=$sample -v mark=$mark -vOFS="\t" '{print sample, mark, $1, $2}' $DMAT >> ${DISTFILE}
    fi
    spref=sampdiff_${sample}_${mark}
    SDMAT=${DIFFDIST_DIR}/${spref}.txt
    if [[ -s $SDMAT ]]; then 
        awk -v sample=$sample -v mark=$mark -vOFS="\t" '{print sample, mark, $1, $2}' $SDMAT >> ${SDISTFILE}
    fi
done < ${TMP_DIR}/imputed_tracks.tsv

gzip -f $DISTFILE
gzip -f $SDISTFILE

# b. Evaluate:



# c. Make tracks for each.
# Lowest distance = nearest.
export RPATH="/broad/compbio/cboix/software/miniconda2/envs/mv_env/bin/R"
source activate mv_env > /dev/null
for mark in `awk 'NR > 1{print $2}' $DISTFILE | sort -u`; do
    echo $mark
    # b. Prioritize strong correlations:
    # NOTE: H3K4me3 comes up a lot. Look at relative to normal.
    # Best way may be a regression - what is the expected similarity of the differential?
    OUTFILE=${TMP_DIR}/prio_${mark}.tsv
    R --slave -f ${BINDIR}/prioritize_diff_tracks.R --args ${DISTFILE} ${mark} ${OUTFILE%%.tsv}

    # Plot best/worst candidates overall by difference:
    MARKFILE=${TMP_DIR}/${mark}_diff_top_diff.tsv
    sort -k4 -n $OUTFILE > ${TMP_DIR}/${mark}_tmp
    head ${TMP_DIR}/${mark}_tmp -n 10 | awk '{print $1"_"$2, $3, NR}' > $MARKFILE
    sort -u $MARKFILE > ${MARKFILE}.tmp
    awk -vOFS="\t" '{print $1"_"$2, $3}' $SAMPLEMARK_TAB | sort -u | join - ${MARKFILE}.tmp | sort -k4 -n | awk -vOFS="\t" '{gsub("_","\t", $1); print $1,$2, $3}' | awk -vOFS="\t" '{print $1,$2,$3,"impute_"$1"_"$2,"difference_"$1"_"$2, "impute_"$1"_"$4}'> $MARKFILE
    qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=1:30:00 -N impobs_$mark -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; $RPATH --slave -f ${BINDIR}/plot_impobs_tracks.R --args ${MARKFILE} chr1 1500000 2000000"

    # Plot best/worst candidates by quantile (per mark):
    for against in `awk 'NR > 1{print $3}' $OUTFILE | sort -u`; do 
        MARKFILE=${TMP_DIR}/${mark}_${against}_diff_top_quant.tsv
        awk -v against=$against '$3 == against' $OUTFILE | sort -k5 -n > ${TMP_DIR}/${mark}_tmp

        head ${TMP_DIR}/${mark}_tmp -n 10 | awk '{print $1"_"$2, $3, NR}' > $MARKFILE
        sort -u $MARKFILE > ${MARKFILE}.tmp
        awk -vOFS="\t" '{print $1"_"$2, $3}' $SAMPLEMARK_TAB | sort -u | join - ${MARKFILE}.tmp | sort -k4 -n | awk -vOFS="\t" '{gsub("_","\t", $1); print $1,$2, $3}' | awk -vOFS="\t" '{print $1,$2,$3,"impute_"$1"_"$2,"difference_"$1"_"$2, "impute_"$1"_"$4}'> $MARKFILE
        qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=1:30:00 -N impobs_${mark}_${against} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; $RPATH --slave -f ${BINDIR}/plot_impobs_tracks.R --args ${MARKFILE} chr1 1500000 2000000"
    done
done

source deactivate > /dev/null

# --------------------------------------------------------
# FINAL: Also submit within sample different mark 
# for all observed tracks - not differences.
# --------------------------------------------------------
CNUM=$( cut -f1 ${SAMPLEMARK_TAB} | sort -u | wc -l | awk '{print $1}' )

qsub -cwd -t 1-$CNUM -P compbio_lab -l h_vmem=8G -l h_rt=0:45:00 -N withinsample_dist -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_withinsample.sh"

# GOAL: Find places where nearest is not same mark.
# NOTE: Look at ..._metrics.R code

rm -rf ${TMP_DIR}
