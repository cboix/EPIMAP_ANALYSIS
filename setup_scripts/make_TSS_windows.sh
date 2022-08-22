#!/bin/bash
# Create TSS windows using the GENCODE v25 GTF annotation of the human genome

# Parameters to make window files
WINS=( CONTROL2_5k 2_5k 5k 50k 100k 150k 200k )
RANGE=( 2500 2500 5000 50000 100000 150000 200000 )
DIST=( 50000 0 0 0 0 0 0 )
for i in `seq 1 7`
do
    WINDOW=${WINS[$i]} 
    WINRANGE=${RANGE[$i]}
    WINDIST=${DIST[$i]}
    echo "$WINDOW $WINRANGE"

    GENEBED=$GTFDIR/Genes.$WINDOW\.bed
    GSORT=${GTFDIR}/Genes.${WINDOW}_sorted.bed
    if [[ ! -s $GENEBED ]] 
    then
        awk -v win=$WINRANGE -v dist=$WINDIST 'BEGIN{FS="\t"; OFS="\t"}{a = sprintf("%s",$1); b = $4 - win-dist; start=(b > 0?b:0); c=$4+win-dist; end=(c > 0?c:0); print a,start,end,$9,$7}' ${GENE_GTF} > $GENEBED

        if [[ "$WINDOW" == "CONTROL_2_5k" ]] 
        then
            # NOTE So as not to create conflicts, i have edited the lines of Genes.CONTROL2_5k with 0, 0 to read 0 5000 (very few cases)
            awk 'BEGIN{FS="\t"; OFS="\t"}{a =($3 == 0?5000:$3); print $1,$2,a,$4,$5}' $GENEBED > $GENEBED\.tmp
            mv $GENEBED\.tmp $GENEBED
        fi
    fi

    if [[ ! -s $GSORT ]] 
    then
        awk '$1~/chr[0-9XY]+/' ${GENEBED} | sort -k1,1 -k2,2n - > $GSORT
    fi
done

# Make TSS bed file: 
TSSFILE=${GTFDIR}/TSS_of_genes.${GTF_SUFFIX}.tsv
if [[ ! -s $TSSFILE ]] 
then
    awk 'BEGIN{FS="\t"; OFS="\t"}{a = sprintf("%s",$1); gsub("gene_id \"","",$0); gsub("\"; gene_name","\t",$0); print a,$4,$9}' ${GENE_GTF} > $TSSFILE
fi
