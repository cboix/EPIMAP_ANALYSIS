#!/bin/bash
# ====================================================
# Process a full set of project files using snakemake:
# Requires PROJECT and EXPT (CRs ChIP-seq for example)
# ====================================================
export PROJECT=$1  # Project name (for directory)
export EXPT=$2  # Experiment type

# Directories:
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]
then
    export DBDIR=/broad/compbio/cboix/${PROJECT}/db
else
    export DBDIR=$HOME/${PROJECT}/db
fi
export SFTDIR=/broad/compbio/cboix/software/bin
export BINDIR=$HOME/${PROJECT}/bin
export ANNDIR=$DBDIR/Annotation
export ENCDIR=${HOME}/data/EPIMAP_ANALYSIS
export SNDIR=${ENCDIR}/bin/snake_pipe
export MAPPING=${ENCDIR}/db/Annotation/all_submitted_released_biosample_mapping.tsv
mkdir -p $DBDIR/out $ANNDIR $BINDIR
cd $SNDIR

# ---------------------------------------------------
# To standardize, FILEINFO has the following columns:
# link id cell_type epitope assembly output_type
# ---------------------------------------------------
export FILEINFO=${ANNDIR}/${PROJECT}_${EXPT}_bam.tsv
echo "Processing project: ${PROJECT} - ${EXPT} has $( wc -l $FILEINFO | awk '{print $1 - 1}' ) files"

# Run in parallel (for each cell type + epitope combination):
COMBNUM=$( awk -vFS="\t" -vOFS="\t" 'NR>1{print $3,$4}' $FILEINFO | sort -u | wc -l - | awk '{print $1}')
SCHEDOPT="-t 1-$COMBNUM -N ${EPITOPE}_snake -o ${DBDIR}/out/ChIP -hold_jid WCE_snake"

# Command:
qsubcmd="qsub $SCHEDOPT ${BINDIR}/run_snakemake.sh project.snake all"

echo "[STATUS] Running the qsub array job:" 
echo "$qsubcmd"
bash -c "$qsubcmd"
