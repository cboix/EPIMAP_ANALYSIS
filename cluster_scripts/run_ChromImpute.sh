#!/bin/bash
# ====================================
# Pass correct sample/mark combination
# to specific ChromImpute command:
# ====================================
if [[ $# -lt 1 ]]; then
    echo "USAGE: $(basename $0) [COMMAND] (optional [INFOFILE] [TASK])" >&2
    echo '  [COMMAND]: ChromImpute command to execute out of:' >&2
    echo '             DISTANCE, FEATURES, PREDICTORS, IMPUTE, EVAL, EXPORT, QEXPORT' >&2
    echo '  [TASK]: (OPTIONAL) sample/mark combination to run (line number from info table)' >&2
    exit 1
fi
COMMAND=$1
start=`date +%s`
hostname -f

if [[ $# -gt 1 ]]; then
    INFOFILE=$2
else
    INFOFILE=${FULLIMPUTATION_TAB}
fi

if [[ $# -gt 2 ]]; then
    SAMPLETABLE=$3
else
    SAMPLETABLE=${SAMPLEMARK_TAB}
fi

if [[ $# -gt 3 ]]; then
    TASK=$4
else
    TASK=${SGE_TASK_ID}
fi

TMP_DIR=${TMP}/${COMMAND}_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR}

if [[ "$CHROM" != "" ]]; then
    CHROMARGS="-c $CHROM"
else
    CHROMARGS=""
fi

echo "USING ARGS: $CHROMARGS"


# Get sample/mark from table (available or to_impute):
if [[ "$COMMAND" == "FEATURES" ]]; then
    MARK=$( sed "${TASK}q;d" ${MARKS_LIST} | awk -v FS="\t" '{print $1}' )
    echo "ChromImpute: Running $COMMAND on $MARK"
elif [[ "$COMMAND" == "DISTANCE" ]] || [[ "$COMMAND" == "EVAL" ]]; then
    SAMPLE=$( sed "${TASK}q;d" ${SAMPLETABLE} | awk -v FS="\t" '{print $1}' )
    MARK=$( sed "${TASK}q;d" ${SAMPLETABLE} | awk -v FS="\t" '{print $2}' )
    INPUTFILE=$( sed "${TASK}q;d" ${SAMPLETABLE} | awk -v FS="\t" '{print $3}' )
    echo "ChromImpute: Running $COMMAND on $SAMPLE and $MARK"
else
    SAMPLE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
    MARK=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $2}' )
    echo "ChromImpute: Running $COMMAND on $SAMPLE and $MARK"
fi

# ChromImpute commands:
if [[ "$COMMAND" == "DISTANCE" ]]; then 

    if [[ ! -s ${DISTANCE_DIR}/${SAMPLE}_${MARK}.txt ]]; then
        echo "2. Compute global distance between datasets with ComputeGlobalDist"
        java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} ComputeGlobalDist -s ${SAMPLE} ${MARK} ${CONVERTED_DIR} ${SAMPLETABLE} ${CHROMSIZES} ${DISTANCE_DIR}
    else
        echo "2. Global distance between datasets already exists. Exiting."
    fi

elif [[ "$COMMAND" == "FEATURES" ]]; then 

    echo "3. Generate training features with GenerateTrainData"
    java -mx35000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} GenerateTrainData ${CHROMARGS} ${CONVERTED_DIR} ${DISTANCE_DIR} ${SAMPLETABLE} ${CHROMSIZES} ${TRAIN_DIR} ${MARK}

elif [[ "$COMMAND" == "PREDICTORS" ]]; then 

    echo "4. Generate the trained predictors for a specific mark + sample type with Train"
    java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} Train ${TRAIN_DIR} ${SAMPLETABLE} ${PREDICTOR_DIR} ${SAMPLE} ${MARK}

elif [[ "$COMMAND" == "IMPUTE" ]]; then 

    cd $PREDICTOR_DIR
    IMPUTED_SUFFIX=impute_${SAMPLE}_${MARK}.wig.gz

    # Remove all old files ONLY if running all at a time.
    if [[ "$CHROM" == "" ]]; then
        echo "TASK $TASK - checking date on files for $SAMPLE and $MARK"
        # TODO: Generalize - shouldn't just be #### May #### 
        rmcmd="rm $( ls -clt $( grep ${SAMPLE}_${MARK} ${CIDIR}/predictor_dir_files ) | awk '$0 !~/ May / && $0 !~ / July / && $0 !~ / August /{print $9}' | awk '{printf $1" "}' )"
        if [[ "$rmcmd" != "rm " ]];then
            echo "PWD: $( pwd )"
            echo "Removing files..."
            echo "COMMAND IS: $rmcmd"
            bash -c "$rmcmd"
        fi
        # Check all chromosomes:
        for file in `ls ${IMPUTED_DIR}/chr*_${IMPUTED_SUFFIX}`; do
            gzip -t $file
            EXITCODE=$?
            if [[ "$EXITCODE" != "0" ]]; then
                echo "Removing $file"
                rm $file
            fi
        done
    else
        echo "TASK $TASK - checking completion of file for $SAMPLE and $MARK - $CHROM:"
        file=${IMPUTED_DIR}/${CHROM}_${IMPUTED_SUFFIX}
        gzip -t $file
        EXITCODE=$?
        if [[ "$EXITCODE" != "0" ]]; then
            echo "Removing $file"
            rm $file
        fi
    fi

    # NOTE: Requires memory to process faster. (othw takes ~24hrs per imputation)
    # Can take NUMLINES up to 100,000 for mx48000M, but that is actually slower.
    echo "5. Generate the imputed signal track for a mark + sample with Apply"
    java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} Apply ${CHROMARGS} ${CONVERTED_DIR} ${DISTANCE_DIR} ${PREDICTOR_DIR} ${SAMPLETABLE} ${CHROMSIZES} ${IMPUTED_DIR} ${SAMPLE} ${MARK}

elif [[ "$COMMAND" == "COEFFV" ]]; then 
    cd $PREDICTOR_DIR
    IMPUTED_SUFFIX=impute_${SAMPLE}_${MARK}_coeffv.wig.gz

    # Remove all old files ONLY if running all at a time.
    echo "TASK $TASK - checking completion of file for $SAMPLE and $MARK - $CHROM:"
    COEFFV_DIR=$CIDIR/coeffv
    mkdir -p $COEFFV_DIR
    file=${COEFFV_DIR}/${CHROM}_${IMPUTED_SUFFIX}
    gzip -t $file
    EXITCODE=$?
    if [[ "$EXITCODE" != "0" ]]; then
        echo "Removing $file"
        rm $file
    fi

    # NOTE: Requires memory to process faster. (othw takes ~24hrs per imputation)
    # Can take NUMLINES up to 100,000 for mx48000M, but that is actually slower.
    echo "5b. Generate the coefficient of variation signal track (-coeffv) for a mark + sample with Apply"
    cmd="java -mx16000M -Djava.io.tmpdir=${TMP_DIR} -jar $BINDIR/ChromImpute_SRC/ChromImpute_coeffv.jar Apply ${CHROMARGS} -coeffv ${CONVERTED_DIR} ${DISTANCE_DIR} ${PREDICTOR_DIR} ${SAMPLETABLE} ${CHROMSIZES} ${COEFFV_DIR} ${SAMPLE} ${MARK}"
    echo "$cmd"
    bash -c "$cmd"

elif [[ "$COMMAND" == "IMPUTE_crossSample" ]]; then 
    # Impute without using cross 
    cd $PREDICTOR_DIR
    IMPUTED_SUFFIX=impute_${SAMPLE}_${MARK}.wig.gz

    echo "5b. Generate the imputed signal track for a mark + sample with Apply without the within-sample information"
    java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} Apply ${CHROMARGS} ${CONVERTED_DIR} ${DISTANCE_DIR} ${PREDICTOR_DIR} ${SAMPLETABLE} ${CHROMSIZES} ${IMPSAMPLE_DIR} ${SAMPLE} ${MARK}
elif [[ "$COMMAND" == "EVAL" ]]; then
    CONVERTED_FILE=${INPUTFILE}.wig.gz
    IMPUTED_FILE=impute_${SAMPLE}_${MARK}.wig.gz
    OUTPUT_FILE=${EVAL_DIR}/${SAMPLE}_${MARK}_eval.txt

    echo "6. Comparing imputed signal track against the converted track with Eval"
    java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} Eval -o ${OUTPUT_FILE} ${CONVERTED_DIR} ${CONVERTED_FILE} ${IMPUTED_DIR} ${IMPUTED_FILE} ${CHROMSIZES}

elif [[ "$COMMAND" == "EXPORT" ]]; then
    inputinfofile=${TMP}/${SAMPLE}_${MARK}_infoimp_${RANDOM}.tsv
    IMPUTED_FILE=impute_${SAMPLE}_${MARK}.wig.gz
    echo "${SAMPLE}	${MARK}	${IMPUTED_FILE}" > ${inputinfofile}
    cat ${inputinfofile}

    mkdir -p ${IMPFMT_DIR}/${MARK}

    # NOTE: Binarize with signalthresh = 2 as stated in paper. Takes max of 10 min.
    echo "7. Convert the imputed signal track for a mark + sample with ExportToChromHMM"
    java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} ExportToChromHMM -b "200" -g "2" "${IMPUTED_DIR}" "${inputinfofile}" "${CHROMSIZES}" "${IMPFMT_DIR}/${MARK}"

    rm ${inputinfofile}
elif [[ "$COMMAND" == "QEXPORT" ]]; then
    # Choose appropriate cutoff:
    CUTFILE=${CIDIR}/binarization_cutoffs_match_avg_distr.tsv
    CUTOFF=$( awk -v mark=$MARK '$1 == mark{print $2}' $CUTFILE )
    echo "[STATUS] Using cutoff $CUTOFF"
    inputinfofile=${TMP}/${SAMPLE}_${MARK}_infoimp_${RANDOM}.tsv
    IMPUTED_FILE=impute_${SAMPLE}_${MARK}.wig.gz
    echo "${SAMPLE}	${MARK}	${IMPUTED_FILE}" > ${inputinfofile}
    cat ${inputinfofile}

    mkdir -p ${IMPFMTQCUT_DIR}/${MARK}

    # NOTE: Binarize with signalthresh = 2 as stated in paper. Takes max of 10 min.
    echo "7 (alt). Convert the imputed signal track for a mark + sample with ExportToChromHMM"
    java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${CHROMIMPUTE} ExportToChromHMM -b "200" -g "$CUTOFF" "${IMPUTED_DIR}" "${inputinfofile}" "${CHROMSIZES}" "${IMPFMTQCUT_DIR}/${MARK}"

    rm ${inputinfofile}
else
    echo "Unknown ChromImpute command: $COMMAND"
    echo "Please use one of: DISTANCE, FEATURES, PREDICTORS, IMPUTE, EVAL, EXPORT"
    exit
fi

rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished $COMMAND sucessfully in $runtime seconds."
