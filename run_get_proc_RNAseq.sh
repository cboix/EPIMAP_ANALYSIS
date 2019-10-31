#!/bin/bash
# --------------------------------------------
# Get specific RNA-seq tsv dataset and process
# --------------------------------------------
start=`date +%s`
hostname -f

# Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES

export TASK=${SGE_TASK_ID}
echo "[STATUS] TASK NUMBER $TASK"

# Directories:
export RNADIR=$DBDIR/RNA-seq
export RIDIR=$RNADIR/file_links/all_submitted_released/
export RNAINFO=$RIDIR/RNA-seq_tsv.csv
export FDIR=$RNADIR/files/RNA-seq
mkdir -p $RNADIR/files $FDIR $FDIR/tsv $FDIR/qc

# Get relevant info:
RURL=$( sed "${TASK}q;d" ${RNAINFO} | awk -F"\t" '{print $1}' )
FILENAME=$( sed "${TASK}q;d" ${RNAINFO} | awk -F"\t" '{sub(".*/","",$1); print $1}' )
ASSEMBLY=$( sed "${TASK}q;d" ${RNAINFO} | awk -F"\t" '{print $3}' )
TYPE=$( sed "${TASK}q;d" ${RNAINFO} | awk -F"\t" '{print $4}' )
BSSID=$( sed "${TASK}q;d" ${RNAINFO} | awk -F"\t" '{print $5}' )
CELL=$( sed "${TASK}q;d" ${RNAINFO} | awk -F"\t" '{print $7}' )

# Files:
ftag=${BSSID}_${FILENAME%%.tsv}
OUTFILE=$FDIR/tsv/${ftag}.tsv
FPKMFILE=$FDIR/tsv/${ftag}_fpkm.tsv

echo "[STATUS] Getting file: $FILENAME ($ASSEMBLY) for id $BSSID - $CELL"
if [[ ! -s ${FPKMFILE} ]]; then
    if [[ ! -s ${OUTFILE}.gz ]]; then
        echo "[STATUS] Downloading file"
        wget $RURL -O $OUTFILE --quiet
        gzip -f $OUTFILE
    else
        echo "[STATUS] File already exists"
    fi
    # Extract names, length, FPKM:
    if [[ "$TYPE" == "transcript quantifications" ]]; then
        echo "[STATUS] Processing $TYPE as transcripts file"
        zcat ${OUTFILE}.gz | awk -vOFS="\t" -F"\t" '$2 ~ /ENSG/{print $2,$7,$3,$1}' | sort +0 -1 > $FPKMFILE
    else
        echo "[STATUS] Processing $TYPE as gene file"
        zcat ${OUTFILE}.gz | awk -vOFS="\t" -F"\t" '$1 !~/^#/ && $1 ~ /ENSG/{print $1,$7}' | sort +0 -1 > $FPKMFILE
    fi
else
    echo "[STATUS] Processed file already exists"
fi

end=`date +%s`
runtime=$((end-start))
echo "Finished RNA-seq download and process sucessfully in $runtime seconds."
