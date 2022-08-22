#!/bin/bash
# ----------------------------------
# Intersect each gwas with an index
# of regions, and report top SNP per
# region + create matrix
# ----------------------------------
GWLIST=$1
SORTIND=$2

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 15

# Make directory:
TMP_DIR=${TMP}/intersect_gid_${SGE_TASK_ID}_${RANDOM}
mkdir -p ${TMP_DIR}

# Get file + attr:
FILE=$( awk -v line=$SGE_TASK_ID 'NR==line{print $1}' $GWLIST )
FNAME=$(basename $FILE)
GWNAME=${FILE##*/gwas/}
GWNAME=${GWNAME%%/*}
BASE=${FNAME%%.bed.gz}
GWINTDIR=$DML_DIR/gwas/$GWNAME
mkdir -p $GWINTDIR

# Run intersection:
echo "[STATUS] Intersection in $GWNAME for trait $BASE"

# NOTE: Pull appropriate p or pval column
zcat ${FILE} | awk -vOFS="\t" 'NR==1{tag=($0~/pval/?"pval":"p"); for(i=1;i<=NF;i++){ix[$i]=i}}NR>1{print "chr"$1, $2, $3 ,$ix[tag]}' > ${TMP_DIR}/tmp_file.tsv

head ${TMP_DIR}/tmp_file.tsv | awk 'NR==1{print "Extracted "NF" fields"}' 

bedtools intersect -loj -sorted -a $SORTIND -b ${TMP_DIR}/tmp_file.tsv > ${TMP_DIR}/tmp_${BASE}.bed

OUTFILE=${GWINTDIR}/$BASE.intersect.tsv.gz
OUTTWO=${GWINTDIR}/$BASE.t1e2.intersect.tsv.gz

awk -vOFS="\t" '$7 != -1{print $4, $8}' ${TMP_DIR}/tmp_${BASE}.bed | sort | awk -vOFS="\t" 'NR==1{n=$1;p=$2}NR > 1{if($1 != n){print n,p; n=$1;p=$2} else {a=(p>$1?$1:p); p=a}}END{print n,p}' | gzip -c > $OUTFILE

zcat $OUTFILE | awk '$2 < 1e-2' | gzip -c > $OUTTWO

# Query lines and filesize
zcat $OUTFILE | wc -l
zcat $OUTTWO | wc -l
ls -sh ${OUTFILE}

rm -rf $TMP_DIR
