#!/bin/bash
# ----------------------------------------
# GWAS hypergeometric enrichments pipeline
# ----------------------------------------
# BEDFILE="/broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust/cls_merge2_wH3K27ac100_nonovl_300_seed1_assignments.fixed.bed"
# RESULTDIR="/broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust/cls_merge2_wH3K27ac100_nonovl_300_seed1/"

# For epigenomes:
# BEDFILE=/broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust/cls_merge2_wH3K27ac100_raw/motif_enrichment/observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_full.bed.gz
# RESULTDIR=/broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust/cls_merge2_wH3K27ac100_raw/

# 0. Load in variables:
BEDFILE=""
RESULTDIR=""
UQFILE=""
FILEPREFIX=""
TAGLINE=""
if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $(basename $0) -i [INPUT BEDFILE] -o [RESULTDIR] (OTHERARGS)
    -i  Input bedfile [REQUIRED] 
            Must have locations + name of set/cluster
    -o  Output directory for results [REQUIRED]
    -u  Unique regions file (e.g. coords) mapping with input (OPTIONAL)
    -f  Fileprefix for centers or names (for plotting, use with CLUSTER/EXTRACT)
    -t  Title/tag (OPTIONAL)"
    exit 1
fi

while getopts i:o:u:f:t: o
do      case "$o" in
    i)		BEDFILE="$OPTARG";;
    o)      RESULTDIR="$OPTARG";;
    u)      UQFILE="$OPTARG";;
    f)      FILEPREFIX="$OPTARG";;
    t)      TAGLINE="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

echo "[STATUS] Parsing opts"

if [[ "$BEDFILE" == "" ]] || [[ "$RESULTDIR" == "" ]]; then
    echo "Did not provide a bedfile or an output directory"
    exit 1
fi

echo "[STATUS] BEDFILE: $BEDFILE"
echo "[STATUS] OUTPUT DIRECTORY: $RESULTDIR"

# NOTE: Load after case catches variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
start=`date +%s` # For timing the pipeline

# Directories 
TMP_DIR=$TMP/run_gwas_analysis_${RANDOM}
mkdir -p ${TMP_DIR} 

# -------------------------------
# 1. Run the enrichments analysis
# -------------------------------
PROCFILE=${RESULTDIR}/hg_processed_enr.Rda
PROCTSV=${RESULTDIR}/hg_processed_enr_long.tsv.gz

# NOTE: need to set path for Rscript and R to work (conflicting R versions)
export RPATH="/broad/compbio/cboix/software/miniconda3/envs/mv_env/bin/"

# conda activate mv_env
if [[ ! -s $PROCFILE ]] || [[ ! -s $PROCTSV ]]; then
    echo "[STATUS] Calculating enrichments"
    # cmd="R --slave -f $BINDIR/gwas_enrichment/calculate_flat_hg_enr.R --args $BEDFILE $RESULTDIR $USE_CORECOORD '$TAGLINE'"
    cmd="$RPATH/Rscript $BINDIR/gwas_enrichment/calculate_flat_hg_enr.R --bedfile $BEDFILE --resultdir $RESULTDIR --uqfile $UQFILE --tagline '$TAGLINE'"
    echo "$cmd"
    bash -c "$cmd"
fi 


# --------------------------------
# 2. Plot the enrichments analysis
# --------------------------------
if [[ "$FILEPREFIX" != "" ]]; then
    # TODO: Plots of stats and overall enrichments
    # TODO: Highlight some of the top enrichments, some of the bottom enrichments.
    echo "[STATUS] Plotting gwas enrichment results"
    cmd="$RPATH/R --slave -f $BINDIR/gwas_enrichment/plot_gwas_alone.R --args $PROCTSV $FILEPREFIX '$TAGLINE' 2500 $RESULTDIR"
    echo "$cmd"
    bash -c "$cmd"
fi


rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished submitting gwas analysis pipeline in $runtime seconds."
