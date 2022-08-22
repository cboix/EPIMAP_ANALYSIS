#!/bin/bash
# Put together all expression datasets
EXPDIR=${DBDIR}/expression
EXPFILE=$EXPDIR/expression_all_FPKM.tsv
cd $EXPDIR

awk '{print $1}' ${EXPDIR}/ENCFF004HYK.tsv > ${EXPDIR}/template_genecol # GRCh38 template
cat $EXPDIR/template_genecol > $EXPFILE

while read file
do
    # ALL or ONLY TAKE GRCh38??
    f2=${file##*/}
    fid=${f2%.*}
    echo $fid
    awk '{print $1,$7}' $file | join ${EXPDIR}/template_genecol - -a 1 -e 0 | awk -v fid=$fid 'BEGIN{print fid}NR > 1{print $2}' | pr -mts $EXPFILE - > ${EXPFILE}.tmp
    mv ${EXPFILE}.tmp $EXPFILE

done < <( ls ${EXPDIR}/ENCFF*.tsv )
