#!/bin/bash
# -----------------------------------
# Process single file bed -> sumstats
# -----------------------------------
if [[ $# -lt 1 ]]; then
    echo "USAGE: $(basename $0) [INFOFILE] (optional [TASK])" >&2
    echo '  [INFOFILE]: GWAS list (file subdir trait)' >&2
    echo '  [TASK]: (OPTIONAL) line number from info table to run' >&2
    exit 1
fi
INFOFILE=$1
start=`date +%s`
hostname -f

if [[ $# -gt 1 ]]; then
    TASK=$2
else
    TASK=${SGE_TASK_ID}
fi

# General variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
export LDSCDIR=$SFTDIR/ldsc

# Get vars for specific task:
file=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
subdir=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $2}' )
trait=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $3}' )
echo "[STATUS] Processing trait $trait in subdir ${subdir}"


TMP_DIR=${TMP}/process_gwas_sumstats_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR}

mkdir -p ${MAINGWAS}/${subdir}/sumstats/
SUMPREF=${MAINGWAS}/${subdir}/sumstats/${trait}
SUMFILE=${SUMPREF}.sumstats.gz
# Check file:
if [[ -s $SUMFILE ]]; then
    gzip -t $SUMFILE
    GZCHECK=$?
    if [[ "$GZCHECK" != "0" ]]; then
        rm $SUMFILE
    fi
fi

# Turn into sumstats:
if [[ ! -s $SUMFILE ]]; then
    NEST=50000
    # TODO determine appropriate # for N, keep consistent throughout
    source activate ldsc
    python $LDSCDIR/munge_sumstats.py --sumstats ${file}\
        --merge-alleles $H3SNPS \
        --out ${SUMPREF} \
        --N $NEST # --a1-inc # In example but not common?
    source deactivate
fi

rm -rf ${TMP_DIR}
end=`date +%s`
runtime=$((end-start))
echo "Finished processing in $runtime seconds."
