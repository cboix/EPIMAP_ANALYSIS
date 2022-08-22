#!/bin/bash
# ---------------------------------------------
# Run ChromHMM on experimental and imputed data
# ChromHMM code by J.Ernst:
# http://compbio.mit.edu/ChromHMM/ 
# ---------------------------------------------
NUMSTATES=18
MFILE=""
MNAME=""
MPREF=""
QCUT="0"
if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [OPTIONS] 
    -i     Info table [required]
    -n     Number of states model (15, 18, 25) (OPTIONAL) [default: $NUMSTATES]
    -p     Prefix to use - to distinguish data used (OPTIONAL)
    -c     Use specific chromHMM file (OPTIONAL) [defaults are the Roadmap models]
    -q     Use matched cutoffs (OPTIONAL - 0/1) [default: $QCUT]
    -m     Model name (prefix) for specific chromHMM file (OPTIONAL, required if using -c)"
    exit 1
fi

while getopts i:n:p:c:m:q: o
do      case "$o" in
    i)		INFOFILE="$OPTARG";;
    n)		NUMSTATES="$OPTARG";;
    p)		MPREF="$OPTARG";;
    c)		MFILE="$OPTARG";;
    q)		QCUT="$OPTARG";;
    m)      MNAME="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES

# Change model, if specifically given one:
# NOTE: epitopes are set as NSTATES right now.

if [[ "$MFILE" != "" ]] && [[ "$MNAME" != "" ]]; then
    export MODEL_LOC=${MFILE}
    export MODELNAME=${MNAME}
fi

if [[ "$MPREF" != "" ]];then
    export MODELPREF="${MODELNAME}_on_$MPREF"
else 
    export MODELPREF=${MODELNAME}
fi

if [[ "$QCUT" == "1" ]]; then
    export MODELPREF="${MODELPREF}_QCUT"
    export FMTDATA_DIR=${IMPFMTQCUT_DIR}
else
    export FMTDATA_DIR=${IMPFMT_DIR}
fi

# Set model-specific directories:
export NCALLDIR=${CALLDIR}/${MODELPREF}
export STATEDIR=${NCALLDIR}/STATEBYLINE
export MODEL_FILE=${ANNDIR}/ChromHMM_model_${MODELNAME}_states.txt
mkdir -p ${NCALLDIR} ${STATEDIR} ${FREQDIR}
cat $MODEL_LOC > $MODEL_FILE # TODO: CLEAN UP OPT.

# Cell from cellinfo table:
cell=$( cut -f1 $INFOFILE | sort -u | sed "${SGE_TASK_ID}q;d" )
CELL_DATA_DIR=${BYCELLDIR}/${cell} 
echo "[STATUS] Running ChromHMM model ${MODELPREF} on sample $( grep ${cell} ${MAPFILE} | awk '{print $2": "$1}' )"
echo "[STATUS] Using $INFOFILE as info file"
mkdir -p ${CELL_DATA_DIR}

# Files/checkpoints:
export FILELIST="${NCALLDIR}/${cell}_${MODELPREF}states_filelist.txt";
export OUTCALLS=${NCALLDIR}/${cell}_${NUMSTATES}_CALLS_segments.bed
export STFILE=${STATEDIR}/${cell}_${NUMSTATES}_CALLS_PER_LINE_chr1_statebyline.txt.gz 
export FREQTAB=${FREQDIR}/${MODELPREF}/table_${cell}.txt

# Ensure frequency table is ok:
if [[ -s $FREQTAB ]]; then
    FTCHECK=$( awk -v ns=$NUMSTATES '$1==ns{print $2 > 100000}' $FREQTAB )
    if [[ "$FTCHECK" != "1" ]]; then
        rm $FREQTAB
    fi
fi

if [[ ! -s $FREQTAB ]]; then
    # ----------------------------------
    # Copy experimental and imputed data
    # to create a complete set of marks
    # ----------------------------------
    # NOTE: Clean all data for the cell on each run 
    # This may interfere with running different models at the same time.
    rm ${CELL_DATA_DIR}/${cell}_${MODELPREF}_chr*binary.txt*
    rm ${CELL_DATA_DIR}/${cell}_${MODELPREF}_origin.txt
    while read printep; do
        if [[ "$printep" == "DNase" ]]; then
            epitope="DNase-seq"
        elif [[ "$printep" == "H2A.Z" ]]; then
            epitope="H2AFZ"
        else 
            epitope=$printep
        fi
        suffix=$( grep -P "${cell}\t${epitope}\t" $INFOFILE | awk '{print $3}')
        echo "$suffix"
        if [[ "$suffix" != "" ]];then
            while read chr size; do
                SUFFIX=${chr}_binary.txt
                EXPFILE=${CHMM_FMTDIR}/${epitope}/${cell}_${DATATAG}_${SUFFIX}.gz
                IMPFILE=${FMTDATA_DIR}/${epitope}/${cell}_${SUFFIX}.gz
                OUTFILE=${CELL_DATA_DIR}/${cell}_${MODELPREF}_${SUFFIX}
                # Add indicated track:
                if [[ "$suffix" == "impute_${cell}_${epitope}" ]]; then
                    usefile=${IMPFILE}
                    usetype="imputed"
                else 
                    usefile=${EXPFILE}
                    usetype="experiment"
                fi
                # NOTE: Files must be column bound together.
                if [[ ! -s ${usefile} ]]; then
                    echo "[STATUS] Indicated ($usetype) file DOES NOT EXIST for ${cell} ${epitope} in ${chr}. EXITING."
                    exit 1
                fi
                echo "[STATUS] Adding ${usetype} file: ${usefile}"
                echo "${printep}\t${chr}\t${usetype}" >> ${CELL_DATA_DIR}/${cell}_${MODELPREF}_origin.txt
                if [[ -s $OUTFILE ]];then
                    zcat $usefile | awk 'NR==1{print ""}NR>1{print $0}' | pr -mts $OUTFILE - > ${OUTFILE}\_tmp
                    mv ${OUTFILE}\_tmp $OUTFILE
                else 
                    zcat $usefile | awk -v head="${cell}\t${chr}" 'NR==1{print head}NR>1{print $0}' > $OUTFILE 
                fi
            done < ${CHROMSIZES_noY}
        else
            echo "[STATUS]"
        fi
    done < <(awk -F"\t" '$1 == "emissionprobs" && $2=="1" && $5 == "0"{print $4}' $MODEL_LOC)


    # Compress all and make filelist:
    gzip -f ${CELL_DATA_DIR}/${cell}_${MODELPREF}_chr*_binary.txt
    find ${CELL_DATA_DIR} -name "${cell}_${MODELPREF}_chr*_binary.txt.gz" -exec basename {} \; | sort > ${FILELIST}

    # -------------------------------------
    # Call chromatin states on binary data:
    # -------------------------------------
    if [[ ! -s ${OUTCALLS}.gz ]] || [[ "$( LC_ALL=C gzip -l ${OUTCALLS}.gz | awk 'NR==2{print $2}')" == "0" ]]; then
        echo "Calling states to ${OUTCALLS}."
        java -mx$((MEM * 1024))M -jar ${CHROMHMM} MakeSegmentation -b 200 -f ${FILELIST} -i CALLS -l ${CHROMSIZES_noY} ${MODEL_FILE} ${CELL_DATA_DIR} ${NCALLDIR}

        # Check file and compress it:
        echo "Checking and compressing output of ${OUTCALLS}."
        stat $OUTCALLS
        sort -k1,1V -k2,2n ${OUTCALLS} | gzip -c > ${OUTCALLS}.gz
    fi

    # ------------------------------------------------
    # Call the chromatin states in STATEBYLINE format:
    # ------------------------------------------------
    if [[ ! -s ${STFILE} ]] || [[ "$( LC_ALL=C gzip -l ${STFILE} | awk 'NR==2{print $2}')" == "0" ]]; then
        echo "Calling states by line."
        java -mx$((MEM * 1024))M -jar ${CHROMHMM} MakeSegmentation -b 200 -nobed -printstatesbyline -f ${FILELIST} -i CALLS_PER_LINE -l ${CHROMSIZES_noY} ${MODEL_FILE} ${CELL_DATA_DIR} ${NCALLDIR}
        gzip -f ${STATEDIR}/${cell}_${NUMSTATES}_CALLS_PER_LINE_*_statebyline.txt
    fi

    # ----------------------------
    # Calculate model frequencies:
    # ----------------------------
    mkdir -p ${FREQDIR}/${MODELPREF};
    if [[ ! -s $FREQTAB ]] || [[ $STFILE -nt $FREQTAB ]]; then
        echo "[STATUS] Calculating state frequencies."
        zcat ${STATEDIR}/${cell}_${NUMSTATES}_CALLS_PER_LINE_*_statebyline.txt.gz | grep "^[0-9]" | \
            sort -n | uniq -c | awk '{print $2, $1}' > ${FREQTAB}
    fi

    echo "Chromatin state frequencies:";
    cat $FREQTAB 
fi
