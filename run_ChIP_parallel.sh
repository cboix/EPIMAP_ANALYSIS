#!/bin/bash
# Start SGE Array: Process ChIP-seq epitopes.
# Note should use process large batch epitopes only
export EPITOPE=$1
export FILEINFO=${CHPLNK}/${EPITOPE}_bam.csv
export EPDIR=${CHPDIR}/files/${EPITOPE}
export TADIR=${EPDIR}/tagAlign
export PKDIR=${EPDIR}/peaks
export BDGDIR=${EPDIR}/bedgraph
export QCDIR=${EPDIR}/qc
mkdir -p $EPDIR $TADIR $QCDIR $BDGDIR $PKDIR ${CONVERTED_DATADIR}/${EPITOPE}
echo "Processing epitope: $EPITOPE - has $( wc -l $FILEINFO | awk '{print $1 - 1}' ) files"

# Parallel for each cell type
awk -vFS="\t" 'NR>1{print $5}' $FILEINFO | sort -u > ${EPDIR}/cell_types
CTNUM=$(wc -l ${EPDIR}/cell_types | awk '{print $1}')


# If processing WCE first, hold:
if [[ "$EPITOPE" == "WCE" ]] 
then
    # For testing: 
    # qsub -cwd -t 2 -l h_vmem=30G -N process_${EPITOPE}_parallel -j y -b y -V -r y -o $DBDIR/out/ChIP $BINDIR/run_ChIP_epitope_cell.sh 
    qsub -cwd -t 1-$CTNUM -q long -l h_vmem=20G -N process_${EPITOPE}_parallel -j y -b y -V -r y -o $DBDIR/out/ChIP $BINDIR/run_ChIP_epitope_cell.sh 
else
    qsub -cwd -t 1-$CTNUM -q long -l h_vmem=30G -N process_${EPITOPE}_parallel -j y -b y -V -r y -hold_jid process_WCE_parallel -o  $DBDIR/out/ChIP $BINDIR/run_ChIP_epitope_cell.sh 
fi

