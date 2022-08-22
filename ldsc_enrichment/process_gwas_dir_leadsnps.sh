#!/bin/bash
# -------------------------------------
# Obtain the leading SNPs for each gwas
# -------------------------------------
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

# Directory
GWLDIR=$DBDIR/gwas_leadsnps/
mkdir -p $GWLDIR

# Get vars for specific task:
file=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
subdir=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $2}' )
trait=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $3}' )
echo "[STATUS] Processing trait $trait in subdir ${subdir}"

mkdir -p ${GWLDIR}/${subdir}/

SUMPREF=${GWLDIR}/${subdir}/${trait}
SUMFILE=${SUMPREF}.all.tsv.gz
PRFILE=${SUMPREF}.lead.tsv.gz
# Check file:
if [[ -s $SUMFILE ]]; then
    gzip -t $SUMFILE
    GZCHECK=$?
    if [[ "$GZCHECK" != "0" ]]; then
        rm $SUMFILE
    fi
fi

cd $GWLDIR/${subdir}

# Get the SNPs with at least 10-4,-6,-8,-10:
if [[ ! -s $PRFILE ]]; then
    col=9
    if [[ "$subdir" == "ukb2" ]]; then
        col=11
    fi
    zcat $file | awk -v ii=$col '$ii <= 10^-4 {print $0}' | gzip -c > $SUMFILE
    source activate pytorch_env
    python $BINDIR/ldsc_enrichment/prune_leadsnps.py main --snpfile $SUMFILE --outfile $PRFILE --col $col
    source deactivate
fi

end=`date +%s`
runtime=$((end-start))
echo "Finished processing in $runtime seconds."
