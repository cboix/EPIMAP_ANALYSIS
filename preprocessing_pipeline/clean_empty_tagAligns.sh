#!/bin/bash
# ----------------------------------------
# Script to remove empty tagAlign.gz files
# ----------------------------------------
# TODO: Add gzip -t test for integrity; do in a faster manner (in run_snake?)
# Arguments:
export EPITOPE=$1
if [[ $# -lt 2 ]]; then 
    export PROJECT=EPIMAP_ANALYSIS
else
    export PROJECT=$2
fi

# Directories:
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]
then
    export DBDIR=/broad/compbio/cboix/${PROJECT}/db
else
    export DBDIR=$HOME/${PROJECT}/db
fi

# TAGALIGN File directory:
if [[  "$EPITOPE" == "DNase-seq" ]] || [[  "$EPITOPE" == "ATAC-seq" ]]; then
    export FILEDIR=${DBDIR}/${EPITOPE}/files/${EPITOPE}/tagAlign
else
    export FILEDIR=${DBDIR}/ChIP-seq/files/${EPITOPE}/tagAlign
fi
mkdir -p ${FILEDIR}

echo "[STATUS] Removing empty ${EPITOPE} tagAligns from project ${PROJECT}"

pushd $FILEDIR
while read file
do
    filesize=$( gzip -l $file | awk 'NR==2{print $2}' )
    if (( $filesize == 0 ))
    then
        gzip -l $file
        rm $file
    fi
done < <( ls *tagAlign.gz )
popd # Return to previous directory

# BEDGRAPH File directory:
if [[  "$EPITOPE" == "DNase-seq" ]] || [[  "$EPITOPE" == "ATAC-seq" ]]; then
    export FILEDIR=${DBDIR}/${EPITOPE}/files/${EPITOPE}/bedgraph
else
    export FILEDIR=${DBDIR}/ChIP-seq/files/${EPITOPE}/bedgraph
fi
mkdir -p ${FILEDIR}

echo "[STATUS] Removing empty ${EPITOPE} bedgraphs from project ${PROJECT}"

pushd $FILEDIR
while read file
do
    filesize=$( gzip -l $file | awk 'NR==2{print $2}' )
    if (( $filesize == 0 ))
    then
        gzip -l $file
        rm $file
    fi
done < <( ls *bedgraph.gz )
popd # Return to previous directory
