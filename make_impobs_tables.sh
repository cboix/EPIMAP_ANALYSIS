#!/bin/bash
# Make the imputed/observed tables.
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 25

# ---------------------------------------------
# A. Make samplesheet that prioritizes imputed 
# but fills gaps with observed.
# ---------------------------------------------
if [[ ! -s $IMPOBS_TAB ]]; then
    ls ${IMPUTED_DIR}/chr1_*.gz | awk -v OFS="\t" -v cvdir="${IMPUTED_DIR}/" '$0 ~ /_impute_/{a = $0; sub(cvdir, "", a); sub("chr[0-9XY]*_impute_","",a); 
    sub("\.wig\.gz$","",a);
    print a,"impute_"a}' | sort -u > ${IMPOBS_TAB}

    # Find differences with full impute:
    awk -vOFS="\t" '{print $1"_"$2, $3}' $SAMPLEMARK_TAB | sort -u > ${SAMPLEMARK_TAB}.tmp
    cat $EXTRAIMPUTATION_TAB $FULLIMPUTATION_TAB > $EXTRAIMPUTATION_TAB.tmp
    awk '{print $1"_"$2}' $EXTRAIMPUTATION_TAB.tmp | sort -u | join -v 2 $IMPOBS_TAB - | join ${SAMPLEMARK_TAB}.tmp - | awk -vOFS="\t" '{print $1, $2}' | tee ${CIDIR}/obs_toadd.tsv >> ${IMPOBS_TAB}

    # 2. Copy these specific unimputed tracks to impute dir:
    cd $CONVERTED_DIR
    while read id suffix; do
        if [[ ! -s "${IMPUTED_DIR}/chr10_${suffix}.wig.gz" ]]; then
            echo $id
            cp $( ls chr*_${suffix}* | grep 'chr[0-9XY]*_FINAL_*' ) ${IMPUTED_DIR}
        fi
    done < ${CIDIR}/obs_toadd.tsv

    echo "[STATUS] Copied over observed files to cover missing impute"

    awk -vOFS="\t" '{sub("_","\t", $1); print $0}' $IMPOBS_TAB | sort -k1 -k2  > ${IMPOBS_TAB}.tmp
    mv ${IMPOBS_TAB}.tmp ${IMPOBS_TAB}
    rm ${SAMPLEMARK_TAB}.tmp ${EXTRAIMPUTATION_TAB}.tmp
fi

# ----------------------------------------------
# B. Make samplesheet that prioritizes observed 
# at good QC, but fills gaps with imputed.
# ----------------------------------------------
AGTAB=${DBDIR}/best_worst_impobs_agreement_011119.tsv
if [[ ! -s $MIXOBS_TAB ]]; then
    rm ${MIXOBS_TAB}
    while read cell; do
        for epitope in ${epitopes[@]}; do
            suffix=$( grep -P "${cell}\t${epitope}\t" ${SAMPLEMARK_TAB} | awk '{print $3}' )
            if [[ "${suffix}" != "" ]]; then
                stats=$( grep -P "${cell}\t${epitope}\t" $AGTAB | awk '{print $3}' )
                if [[ "$stats" == "LOW" ]]; then 
                    echo "${suffix} is low."
                    echo "${cell}\t${epitope}\timpute_${cell}_${epitope}" >> ${MIXOBS_TAB}
                else
                    echo "${cell}\t${epitope}\t${suffix}" >> ${MIXOBS_TAB}
                fi
            else
                echo "${cell}\t${epitope}\timpute_${cell}_${epitope}" >> ${MIXOBS_TAB}
            fi
        done
    done < <( cut -f1 ${SAMPLEMARK_TAB} | sort -u )
    sort -k1 -k2 $MIXOBS_TAB > ${MIXOBS_TAB}.tmp
    mv ${MIXOBS_TAB}.tmp ${MIXOBS_TAB}
fi

# -------------------
# All existing tracks
# -------------------
if [[ ! -s ${ALL_TRACKS_TAB} ]] || [[ ! -s ${ALL_UQ_TAB} ]]; then
    # Observed:
    cat $SAMPLEMARK_TAB > ${ALL_TRACKS_TAB}
    while read cell mark; do
        if [[ -s ${IMPUTED_DIR}/chr10_impute_${cell}_${mark}.wig.gz ]]; then
            echo "${cell}\t${mark}\timpute_${cell}_${mark}" >> ${ALL_TRACKS_TAB}
        fi
    done < $FULLIMPUTATION_TAB 

    while read cell mark; do
        if [[ -s ${IMPUTED_DIR}/chr10_impute_${cell}_${mark}.wig.gz ]]; then
            echo "${cell}\t${mark}\timpute_${cell}_${mark}" >> ${ALL_TRACKS_TAB}
        fi
    done < $EXTRAIMPUTATION_TAB
    # With unique cell identifiers for imputed/observed:
    awk -vOFS="\t" '{a = substr($3, 0, 5); suf=(a=="FINAL"?"_obs":"_imp"); print $1suf,$2,$3}' ${ALL_TRACKS_TAB} > ${ALL_UQ_TAB}

    # Copy over all converted files:
    cd $CONVERTED_DIR
    while read cell mark suffix; do
        if [[ ! -s "${IMPUTED_DIR}/chr10_${suffix}.wig.gz" ]]; then
            echo $cell $mark
            cp $( ls chr*_${suffix}* | grep 'chr[0-9XY]*_FINAL_*' ) ${IMPUTED_DIR}
        fi
    done < ${SAMPLEMARK_TAB}
fi
