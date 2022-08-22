#!/bin/bash
# =============================
# Calculates DHS coverage
# Over TSS windows + chrom arms 
# =============================
# Arguments: 
id=$1;
cell=$2;
STRID=${id}_${cell}
TAFILE=${FILE_DIR}/FINAL_${STRID}.tagAlign.gz
TMPFILE=${FILE_DIR}/FINAL_${STRID}.tagAlign.tmp.gz

# Sorted temporary file: 
zcat $TAFILE | sort -k1,1 -k2,2n | gzip -c > $TMPFILE

# Variables:
WINS=( CONTROL2_5k 2_5k 5k 50k 100k 150k 200k )
# TODO MAKE SUBTRACTED WINDOWS -- but after normalizing!

# All windows together per DHS dataset:
DW_FILE=${WINDIR}/${STRID}.DHS_windows
rm -f $DW_FILE
touch $DW_FILE

# Compute window overlaps:
for WIN in ${WINS[@]}
do
    WINBED=$GTFDIR/Genes.${WIN}_sorted.bed
    JOBNAME=$STRID\_$WIN
    TABFILE=$WINDIR/$JOBNAME\.tsv

    echo "Calculating: $JOBNAME"
    bedtools intersect -a $WINBED -b $TMPFILE -c -sorted > $TABFILE 

    # Collate:
    echo "Adding $WIN to main file"
    awk -v vari=$WIN 'BEGIN{FS="\t"; print vari}{print $6}' $TABFILE | pr -mts $DW_FILE - > $DW_FILE\_2 
    mv $DW_FILE\_2 $DW_FILE
done

# Add gene cols: ensembl and header 
echo "Adding gene columns"
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7}' $DW_FILE | sed "1d" | pr -mts ${GENE_COLS} - |  sed "1iCHR\tNAME\tCONTROL2_5k\t2_5k\t5k\t50k\t100k\t150k\t200k"  > $DW_FILE\_2
mv $DW_FILE\_2 $DW_FILE

if [[ "9" == "$( awk 'NR==2{print NF}' $DW_FILE )" ]] 
then
    rm -f ${WINDIR}/${STRID}_*.tsv
fi

# Chromosome occupancy:
CHRFILE=${WINDIR}/${STRID}.chrom_arm
if [[ ! -s $CHRFILE ]] 
then
    bedtools intersect -a $CHROMARM -b $TMPFILE -c > $CHRFILE
fi
rm -f $TMPFILE
