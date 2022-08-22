#!/usr/bin/env bash
# Run enhancer-gene linking:
start=`date +%s`
hostname -f
# Arguments:
export SAMPLE=""
export CHROMHMM_STATE=""
export CUTOFF="2.5"
export DISTANCE="0"
export MARK="0"
export ROADMAP="0"

if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [-t COMMAND] [OPTIONS] 
    -i     Sample [required] [biosample ids in form of: BSSXXXXX]
    -s     ChromHMM state (18-state) to run [required] [E7|E8|E9|E10|E11|E15]
    -c     Cutoff to use when reporting links [default: $CUTOFF]
    -d     Distance at which to boost with distance pred. [default: $DISTANCE]
    -m     Whether to add mark information for prediction (0/1) [default: $MARK]
    -r     Whether to use logistic reg. roadmap method (0/1) [default: $ROADMAP]"
    exit 1
fi

while getopts i:s:c:d:m:r: o
do      case "$o" in
    i)		export SAMPLE="$OPTARG";;
    s)		export CHROMHMM_STATE="$OPTARG";;
    c)		export CUTOFF="$OPTARG";;
    d)      export DISTANCE="$OPTARG";;
    m)      export MARK="$OPTARG";;
    r)      export ROADMAP="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

# Load in variables:
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]; then
    source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
else
    source /home/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
fi
SRCDIR="$BINDIR/eglinking/linking-xgboost/linking/";
DATADIR="$DBDIR/linking_data/";

BASENAME=`echo "${SAMPLE}_${CHROMHMM_STATE}_${CUTOFF}" | tr ' ' '_'`;

echo "Method is $METHOD"

if [[ "${METHOD}" == "ROADMAP" ]]; then
    LINKSRC=$SRCDIR/linking_logistic.py 
else 
    LINKSRC=$SRCDIR/linking.py 
fi

if [[ ! -s $BASENAME.txt ]]; then
    echo "[STATUS] Running predictions to get ${BASENAME}.txt"
    cmd="$LINKSRC --sample $SAMPLE --chromhmm-state $CHROMHMM_STATE \
        -o "${BASENAME}.txt" -l "${BASENAME}.log" 2>"${BASENAME}.err" \
        --dhs-chromhmm $DATADIR/dhs_chromhmm.hdf5 \
        --pos-metadata $DATADIR/pairs_ENH_ovl_df.cp.gz \
        --neg-metadata $DATADIR/random_pairs_ENH_ovl_df.cp.gz \
        --cutoff $CUTOFF \
        --cor $DATADIR/H3K27ac_precomputed_corr.hdf5 \
        $DATADIR/H3K4me1_precomputed_corr.hdf5 \
        $DATADIR/H3K4me2_precomputed_corr.hdf5 \
        $DATADIR/H3K4me3_precomputed_corr.hdf5 \
        $DATADIR/H3K9ac_precomputed_corr.hdf5 \
        $DATADIR/DNase-seq_precomputed_corr.hdf5 \
        --rand-cor $DATADIR/H3K27ac_precomputed_random_corr.hdf5 \
        $DATADIR/H3K4me1_precomputed_random_corr.hdf5 \
        $DATADIR/H3K4me2_precomputed_random_corr.hdf5 \
        $DATADIR/H3K4me3_precomputed_random_corr.hdf5 \
        $DATADIR/H3K9ac_precomputed_random_corr.hdf5 \
        $DATADIR/DNase-seq_precomputed_random_corr.hdf5"

    # Only add other flags if not running roadmap
    if [[ "$ROADMAP" != "1" ]]; then
        # Add the distance flag
        if [[ "$DISTANCE" != "0" ]]; then
            cmd="$cmd --distance-cutoff $DISTANCE" 
        fi

        # Mark 
        if [[ "$MARK" == "1" ]]; then
            cmd="$cmd --mark ${DATADIR}/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5 ${DATADIR}/H3K4me1_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5"
        fi
    fi

    echo "$cmd"

    bash -c "$cmd"

else
    echo "[STATUS] ${BASENAME}.txt already exists"
fi

end=`date +%s`
runtime=$((end-start))
echo "[STATUS] Finished in $runtime seconds."
