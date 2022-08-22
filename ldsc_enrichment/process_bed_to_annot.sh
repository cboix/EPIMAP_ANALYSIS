#!/bin/bash
# ---------------------------------------------------------------
# Process BED files into annot.gz files
# NOTE: Processing follows:
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
# ---------------------------------------------------------------
if [[ $# -lt 1 ]]; then
    echo "USAGE: $(basename $0) [INFOFILE] (optional [TASK])" >&2
    echo '  [INFOFILE]: BED file list (file subdir trait)' >&2
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
BEDFILE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
PREFIX=${bedfile%%.bed}
ANNFILE=${PREFIX}.annot.gz

# Check file integrity:
if [[ -s $ANNFILE ]]; then
    gzip -t $ANNFILE
    GZCHECK=$?
    if [[ "$GZCHECK" != "0" ]]; then
        rm $ANNFILE
    fi
fi

echo "[STATUS] Processing file $BEDFILE"
if [[ ! -s $ANNFILE ]]; then
    python $LDSCDIR/make_annot.py \
        --bed-file $BEDFILE \
        --bimfile 1000G.EUR.QC.22.bim \
        --annot-file $ANNFILE
fi

GENOMEDIR=/broad/compbio/cboix/genomes
E1KDIR=${GENOMEDIR}/ldsc_eur1000kg
HM3DIR=${GENOMEDIR}/hapmap3

# Make the l2 LDscore file from the Annotation:
LASTL2FILE=${PREFIX}.22.l2.ldscore.gz 
if [[ ! -s ${LASTL2FILE} ]]; then
    for chr in `seq 1 22`; do 
        python ${LDSCDIR}/ldsc.py --l2 \ 
            --bfile ${E1KDIR}/1000G.EUR.QC.${chr} \ 
            --ld-wind-cm 1 \ 
            --annot ${ANNFILE} \
            --thin-annot --out ${PREFIX} \
            --print-snps ${HM3DIR}/hm.${chr}.snp
        done
fi

# Note: the --thin-annot flag is only included if the annot file does not have the CHR, BP, SNP, and CM columns.) Repeat for the other chromosomes. Make sure you save your output files to the same directory, with the same file prefix, as your annot files, so that you have ${prefix}.${chr}.annot.gz, ${prefix}.${chr}.l2.ldscore, and ${prefix}.${chr}.l2.M_5_50.


end=`date +%s`
runtime=$((end-start))
echo "Finished processing in $runtime seconds."
