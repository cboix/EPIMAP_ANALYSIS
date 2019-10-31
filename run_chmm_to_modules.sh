#!/bin/bash
# ------------------------------------------
# Pipeline to take ChromHMM into modules
# 0. PREPROCESS
# 1. AGGREGATE
#    Aggregate ChromHMM for elements/states
#    - Intersected with masterlist
# 2. MAKEMAT
#    - Use to merge matrices only
# 3. CLUSTER
#    Cluster resulting matrix to get modules
# 
# NOTE: To run on raw data, flag with -a 1
# and run EXTRACT first to make bedfiles
# ------------------------------------------
start=`date +%s`
hostname -f
# Default model is 18-state, mixed obs/imputed:
export NUMSTATES=18
export MIXOBS="1"
export QCUT="1"
# Other arguments:
export NCLUST="300"
export ELEMENT="ENH"
export SPECIFIC="0"
export EXTSIZE="100"
export PREF="" 
# export SMERGE="0"
# export MWITH="0"
export SMERGE="2"
export MWITH="1"
export RAWSIGNAL="0"
export NOCLS="0"

if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [-t COMMAND] [OPTIONS] 
    -t     COMMAND [required] [one of PREPROCESS, AGGREGATE, MAKEMAT, CLUSTER, etc.]
    -s     Intersect cell-type specific or not (OPTIONAL - 0/1) [default: $SPECIFIC]
    -p     Prefix to use - to distinguish data used (OPTIONAL)
    -e     Element to use (ENH, PROM, state #, or mark/assay) [default: $ELEMENT]
    -k     Number of clusters [default: $NCLUST]
    -g     Group/merge states (OPTIONAL - 0/1/2) [default: $SMERGE]
           where 0 is no merge, 1 is according to sub-type, 2 is full merge
    -r     Use raw signal (for marks/assays) (OPTIONAL - 0/1) [default: $RAWSIGNAL]
    -x     Extension size (for raw data - H3K27ac) [default: $EXTSIZE]
    -w     Merge with H3K27ac (OPTIONAL - 0/1) [default: $MWITH]
    -a     NO CLUSTERING - use raw (OPTIONAL - 0/1) [default: $NOCLS]
    ------------ MODEL INFORMATION ------------
    -n     Number of states model (15, 18, 25) (OPTIONAL) [default: $NUMSTATES]
    -q     Use matched cutoffs (OPTIONAL - 0/1) [default: $QCUT]
    -m     Use mixed obs/imputed data (OPTIONAL - 0/1) [default: $MIXOBS]"
    exit 1
fi

while getopts t:s:p:e:k:g:r:x:w:a:n:q:m: o
do      case "$o" in
    t)		export COMMAND="$OPTARG";;
    s)		export SPECIFIC="$OPTARG";;
    p)		export PREF="$OPTARG";;
    e)      export ELEMENT="$OPTARG";;
    k)      export NCLUST="$OPTARG";;
    g)		export SMERGE="$OPTARG";;
    r)      export RAWSIGNAL="$OPTARG";;
    x)      export EXTSIZE="$OPTARG";;
    w)      export MWITH="$OPTARG";;
    a)      export NOCLS="$OPTARG";;
    n)		export NUMSTATES="$OPTARG";;
    q)		export QCUT="$OPTARG";;
    m)      export MIXOBS="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

# Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES

export TASK=${SGE_TASK_ID}
echo "[STATUS] TASK NUMBER $TASK"

# Specify model:
export MODEL=${MODELNAME}
if [[ "$MIXOBS" == "1" ]]; then
    export INFOFILE=$MIXOBS_TAB
    export MODEL=${MODEL}_on_mixed_impobs
    export JOBSUF="on_mixed_impobs"
else 
    export INFOFILE=$IMPOBS_TAB
    export JOBSUF="on_imputed"
fi

if [[ "$QCUT" == "1" ]]; then
    export MODEL=${MODEL}_QCUT
    export FMTDATA_DIR=${IMPFMTQCUT_DIR}
else
    export FMTDATA_DIR=${IMPFMT_DIR}
fi

echo "[STATUS] Model used is ${MODEL} for element type $ELEMENT"
# Allow creation for other elements
# Allow for specific state numbers

export ISMARK=$( grep -e "^$ELEMENT$" $MARKS_LIST | wc -l )
export CHMM_STATEMAP=${ANNDIR}/chmm_element_mapping.txt
if [[ "$ISMARK" == "1" ]]; then
    export STATES=$ELEMENT
    export TAGLINE="$ELEMENT"
elif [[ "$ELEMENT" == "ENH" ]] ||
    [[ "$ELEMENT" == "PROM" ]]; then
    stline=$( awk -v elem=$ELEMENT -v ns=$NUMSTATES '$1 == elem && $2 == ns' $CHMM_STATEMAP )
    STATES=$( echo $stline | awk '{print $3}' )
    MERGED_STATES=$( echo $stline | awk '{print $4}' )
    export STATES="[$STATES]"
    export MERGED_STATES="[$MERGED_STATES]"
    if [[ "$ELEMENT" == "ENH" ]]; then
        export TAGLINE="ChromHMM Enhancers from $NUMSTATES state model"
    elif [[ "$ELEMENT" == "PROM" ]]; then
        export TAGLINE="ChromHMM Promoters from $NUMSTATES state model"
    fi
else
    export STATES="$ELEMENT"
    export TAGLINE="ChromHMM state $ELEMENT from $NUMSTATES state model"
fi

echo "[STATUS] Set states for element $ELEMENT to $STATES"
echo "[STATUS] Set tagline to: $TAGLINE"

# Make temporary directories:
export TMP_DIR=${TMP}/CTOM_${COMMAND}_${MODEL}_${RANDOM}
mkdir -p ${TMP_DIR}

# --------------------------------------------------------
# Arguments for the matrix intersections and aggregations:
# --------------------------------------------------------
# Create aggregation infofile for specific job
export JOBINFO=${TMP_DIR}/job_infofile.tsv
if [[ "$ISMARK" == "1" ]]; then
    echo "[STATUS] $ELEMENT is mark/assay. Running aggregation regardless of model"
    if [[ "$RAWSIGNAL" == "1" ]]; then
        # NOTE: For imputed data, the raw signal 
        # will be conservative for all but H3K27ac
        echo "[STATUS] Running raw signal aggregation"
        export ELEMDIR=${IMPUTED_DIR}
        export DATADIR=${IMPUTED_DIR}/${ELEMENT}
        export FPREFIX=${ELEMENT}_all_bin_${JOBSUF}
        mkdir -p $DATADIR
        # Populate temporary infofile (id, prefix, suffix)
        awk -vOFS="\t" -v mark=$ELEMENT -v obsdir=$CONVERTED_DIR/ -v impdir=$IMPUTED_DIR/ '$2 == mark{if ($3 ~/impute/){print $1,"","_"$3".wig.gz", impdir} else {print $1,"","_"$3".wig.gz", obsdir}}'  $INFOFILE > $JOBINFO
        export SUBARGS="--noverbose --resolution 25 --extend $EXTSIZE"
    else
        echo "[STATUS] Running binary signal aggregation"
        export ELEMDIR=${CHMM_FMTDIR}/${ELEMENT}
        export DATADIR=${ELEMDIR}
        export FPREFIX=${ELEMENT}_all_bin_${JOBSUF}
        # Populate temporary infofile (id, prefix, suffix)
        awk -vOFS="\t" -v mark=$ELEMENT -v obsdir=$CHMM_FMTDIR/$ELEMENT/ -v impdir=$FMTDATA_DIR/$ELEMENT/ '$2 == mark{if ($3 ~/impute/){print $1, $1"_","_binary.txt.gz", impdir} else {print $1, $1"_c_t.sub_","_binary.txt.gz", obsdir}}' $INFOFILE > $JOBINFO
        export SUBARGS="--noverbose --resolution 200 --extend 0"
    fi
else 
    echo "[STATUS] Aggregating ChromHMM states $STATES (${ELEMENT}) for model ${MODEL}"
    export CALLDIR=${CHMMDIR}/calls/${MODEL}
    export ELEMDIR=${CALLDIR}/STATEBYLINE
    export DATADIR=${CALLDIR}/${ELEMENT}
    export FPREFIX=${MODEL}_${ELEMENT}_bin_${JOBSUF}
    mkdir -p ${DATADIR}
    # Populate temporary infofile (id, prefix, suffix)
    awk -vOFS="\t" -v nstates=$NUMSTATES '{print $1, $1"_"nstates"_CALLS_PER_LINE_","_statebyline.txt.gz"}' $ANNDIR/kept_bssid_20190322.txt > $JOBINFO
    export SUBARGS="--states $STATES --chromhmm --noverbose"
fi
export EOUTPREF=${DATADIR}/$FPREFIX
export MAINARGS="--infofile $JOBINFO --dir ${ELEMDIR}/ --out $EOUTPREF $SUBARGS"
if [[ "$SMERGE" == "2" ]]; then 
    MAINARGS="$MAINARGS --mergestates"
fi


# ---------------------------------------------------
# Set up directories for clusters/motifs/enrichments:
# STRUCTURE: 
# Underneath main clsdir/clsprefix:
# - motif enrich
# - great enrich
# - gwas catalog
# - ldak and ldsc analysis
# - directory plots of each of these analyses
# ---------------------------------------------------
export CLSDIR=${DATADIR}/clust
if [[ "$MWITH" == "1" ]]; then 
    CPREF=cls_merge${SMERGE}_wH3K27ac${EXTSIZE}
else
    CPREF=cls_merge${SMERGE}
fi
if [[ "$NOCLS" == "1" ]]; then
    export CPREFIX=${CPREF}_raw
else
    export CPREFIX=${CPREF}_${NCLUST}
fi
export TCDIR=${CLSDIR}/$CPREFIX  # Top directory
export TCIMG=${TCDIR}/img/
# Sub-directories:
export GTDIR=${TCDIR}/go_enrichment
export GTBKDIR=${TCDIR}/go_enrichment_with_bkg
export MTDIR=${TCDIR}/motif_enrichment
export GWDIR=${TCDIR}/gwas_catalog
export LKDIR=${TCDIR}/gwas_ldak
export LCDIR=${TCDIR}/gwas_ldsc
mkdir -p ${DATADIR} ${CLSDIR} ${TCDIR} ${GTDIR} ${MTDIR} ${GWDIR} ${LKDIR} ${LCDIR} ${IMGDIR} ${GTBKDIR}
# Prefixes/files:
export COUTPREF=${CLSDIR}/$CPREFIX
export CLSBED=${COUTPREF}_assignments.bed
export SPLITPREFIX=${MTDIR}/${FPREFIX}  # Prefix for split files (in mtdir)
export GWPREFIX=${GWDIR}/${FPREFIX}


# ------------------
# Run inner command:
# ------------------
source ${BINDIR}/run_chmm_to_modules_innercmd.sh
# CMLIST=(GREAT GWAS MOTIF FIMO)
# for com in ${CMLIST[@]}; do 
#     export COMMAND=$com
#     source ${BINDIR}/run_chmm_to_modules_innercmd.sh
# done

# Collapse any images that were generated:
cd $TCDIR
tar -cvzf $TCDIR/images.tar.gz img
echo "[STATUS] Images at $TCDIR/images.tar.gz"

rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished $COMMAND sucessfully in $runtime seconds."
