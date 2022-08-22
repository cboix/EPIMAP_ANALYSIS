#!/bin/bash
# ---------------------------------------------------------------
# Process BED files into annot.gz files
# NOTE: Processing follows:
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
# ---------------------------------------------------------------
if [[ $# -lt 1 ]]; then
    echo "USAGE: $(basename $0) [INFOFILE] [SSLIST] [OUTDIR] (optional [TASK])" >&2
    echo '  [INFOFILE]: Annotation BED files (file identifier)' >&2
    echo '  [SSLIST]: Sumstats file list (files gwastype trait)' >&2
    echo '  [OUTDIR]: Output directory' >&2
    echo '  [TASK]: (OPTIONAL) line number from SSLIST to run' >&2
    exit 1
fi
# Necessary arguments:
INFOFILE=$1
SSLIST=$2
RESULTDIR=$3
start=`date +%s`
hostname -f

if [[ $# -gt 3 ]]; then
    TASK=$4
else
    TASK=${SGE_TASK_ID}
fi

# Load general variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
export LDSCDIR=$SFTDIR/ldsc

# Get vars for specific task:
SSFILE=$( sed "${TASK}q;d" ${SSLIST} | awk -v FS="\t" '{print $1}' )
SUBDIR=$( sed "${TASK}q;d" ${SSLIST} | awk -v FS="\t" '{print $2}' )
TRAIT=$( sed "${TASK}q;d" ${SSLIST} | awk -v FS="\t" '{print $3}' )

# Outputs: 
WEIGHTPREF=${LDSCBASE}/weights_hm3_no_hla/weights.
BASEPREF=${LDSCBASE}/1000G_EUR_Phase3_baseline/baseline.
mkdir -p ${RESULTDIR}/${SUBDIR}

# for above:
python ${LDSCDIR}/ldsc.py 
	--h2 $SSFILE \
	--ref-ld-chr ${BASEPREF} \ 
	--w-ld-chr ${WEIGHTPREF} \
	--overlap-annot \  # Tells it that files overlap
	--frqfile-chr 1000G.mac5eur.\
	--out ${RESULTDIR}/${SUBDIR}/${TRAIT}_baseline
    # TODO WHERE MAC5EUR FILE?? 
    # TODO SHOULD WE DO CT SPECIFIC: 

python ${LDSCDIR}/ldsc.py \
    --h2-cts UKBB_BMI.sumstats.gz \
    --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
    --out BMI_${cts_name} \
    --ref-ld-chr-cts $cts_name.ldcts \
    --w-ld-chr weights_hm3_no_hla/weights.

# Uses a testing ldcts file like this: 
# CELLTYPE  TESTPREFIX,CONTROLPREFIX
# CELLTYPE  TESTPREFIX,CONTROLPREFIX

# NOTE: need to come up with matched controls?
# - Avg of node above 
# - Permutations of matrix to preserve row/colsums
# - Weighted (by margin prob) sampling of enhancers to same number N


end=`date +%s`
runtime=$((end-start))
echo "Finished processing in $runtime seconds."
