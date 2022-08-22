#!/bin/bash
# ----------------------------------
# Intersect each gwas with an index
# of regions, and report top SNP per
# region + create matrix
# ----------------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

# Temp dir and gwas output dir:
TMP_DIR=${TMP}/intersect_gid_${RANDOM}
mkdir -p ${TMP_DIR} $DML_DIR/gwas/

# List out all GWAS:
GWASDIR=/broad/compbio/data/gwas
GWLIST=${TMP_DIR}/gwas_bedfiles.tsv
ls $GWASDIR/*/tabix/*.bed.gz > $GWLIST

GNUM=$( cat ${GWLIST} | wc -l ) 
echo "[STATUS] Will intersect against $GNUM traits"

# Make sorted file for faster intersect:
INDFILE=${DMLPREF}_hg19.txt
SORTIND=${DMLPREF}_hg19_sorted.bed
if [[ ! -s $SORTIND ]]; then
    sort -k1,1 -k2,2n $INDFILE > $SORTIND
fi

# Run commands:
qsub -cwd -t 1-$GNUM -P compbio_lab -tc 2500 -l h_vmem=5G -l h_rt=0:30:00 -N int_gwas -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/intersect_gwas_with_index.sh $GWLIST $SORTIND"

# Aggregate all GWAS intersections into matrix
TMPFILE=${TMP_DIR}/tmp_fullintersect.tsv
FULLFILE=$DML_DIR/gwas/all_intersect.t1e2.tsv.gz
while read file; do
    # Get file + attr:
    FNAME=$(basename $file)
    BASE=${FNAME%%.bed.gz}
    GWNAME=${file##*/gwas/}
    GWNAME=${GWNAME%%/*}
    echo $GWNAME $BASE
    GWINTDIR=$DML_DIR/gwas/$GWNAME
    FILTFILE=${GWINTDIR}/$BASE.t1e2.intersect.tsv.gz
    zcat $FILTFILE | awk -vOFS="\t" -v nam="${GWNAME}/${BASE}" '{print nam, $1,$2}' >> $TMPFILE
done < $GWLIST
gzip -c $TMPFILE > $FULLFILE


# Match matrix to clusters matrix


