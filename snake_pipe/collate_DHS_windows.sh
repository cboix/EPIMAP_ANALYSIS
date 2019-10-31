#!/bin/bash
# -----------------------------------------------
# Put together windows (mainly for DHS data)
# and run processing: 
# - combat normalization
# - merge replicates
# - create data matrix
# -----------------------------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

# Arguments:
if [[ $# -lt 1 ]]; then 
    export EPITOPE=DNase-seq
else
    export EPITOPE=$1
fi

# File directory:
if [[  "$EPITOPE" == "DNase-seq" ]] 
then
    export EPDIR=${DBDIR}/${EPITOPE}/files/${EPITOPE}/
    export INFOTABLE=${DBDIR}/${EPITOPE}/file_links/all_submitted_released/${EPITOPE}_bam.csv
else
    export EPDIR=${DBDIR}/ChIP-seq/files/${EPITOPE}/
    export INFOTABLE=${DBDIR}/ChIP-seq/file_links/all_submitted_released/${EPITOPE}_bam.csv
fi

# DNase-seq directories:
export FILEDIR=${EPDIR}/tagAlign
export WINDIR=${EPDIR}/window_coverage
export QCDIR=${EPDIR}/qc
mkdir -p $FILEDIR $QCDIR $WINDIR

# Clean empty tagAligns:
${SNDIR}/clean_empty_tagAligns.sh DNase-seq

# Report expected counts:
cd $EPDIR
FNUM=$(wc -l $INFOTABLE | awk '{print $1}')
CNUM=$(cut -f5 $INFOTABLE | sort -u | wc -l | awk '{print $1}')
echo "[STATUS] There should be ${FNUM} files across ${CNUM} cells"

# Get current count for files (cells) and windows (replicates):
SNUM=$( ls ${FILEDIR} | egrep "FINAL_*_BSS[0-9]+.sub.tagAlign.gz" | wc -l | awk '{print $1}' )
TNUM=$( ls ${FILEDIR} | egrep "FINAL_*_BSS[0-9]+.tagAlign.gz" | wc -l | awk '{print $1}' )
echo "[STATUS] There are ${TNUM} full tagAligns and ${SNUM} subsampled tagAligns (cells)"

# Cleaning + making list for merging window files: 
cd $WINDIR
echo "[STATUS] Cleaning out incomplete window files for ${EPITOPE}"
MINF=${EPDIR}/merging_${EPITOPE}_windows.tsv
rm -f $MINF
cd $WINDIR
while read file
do
    if [[ "9" == "$( awk 'NR==2{print NF}' $file )" ]] 
    then
        id=${file%%_*}
        prefix=${file%.DHS_windows}
        cell=${prefix#*_}
        echo "$file\t$id\t$cell" >> $MINF
    else 
        rm $file
    fi
done < <( ls *.DHS_windows )
WNUM=$( ls ${WINDIR}/*.DHS_windows | wc -l | awk '{print $1}' )
echo "[STATUS] There are ${WNUM} window files"

# Merge files:
for ((i = 3; i <= 9; i++))
do
    echo $i
    # SET UP FILE: 
    WIN=$(awk -v ind=$i 'NR==1{print $ind}' ${WINDIR}/ENCFF402SAN_BSS00004.DHS_windows)
    OUTFILE=${WINDIR}/all_${WIN}_raw.tsv
    if [[ ! -s $OUTFILE ]] 
    then
        echo "STARTING to make $OUTFILE"
        IFS=$'\t'
        while read file id cell
        do
            echo "Add $id, $cell"
            awk -v ind=$i -v id=$id -v cell=$cell 'NR==1{print id"_"cell}NR>1{print $ind}' $file | pr -mts $OUTFILE - > ${OUTFILE}.tmp
            mv ${OUTFILE}.tmp ${OUTFILE}
        done < $MINF
        echo "FINISHED making $OUTFILE"
    fi
done


