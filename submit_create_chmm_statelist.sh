#!/bin/bash
# --------------------------------------------
# Create list of regions for a state 
# from a STATEBYLINE directory
# outputs tsv of chmm bin locations.
# --------------------------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh


# Default is running for n25 states:
STDIR='/broad/compbio/jernst/compbio-hp/SIGNALPREDICT6/TIER2_REORDER_MODEL/states25/STATEBYLINE/'
PREFIX=${DBDIR}/state_regionlists/imputed12model

python ${BINDIR}/create_chmm_state_possiblelist.py main --directory ${STDIR} --outprefix $PREFIX

qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=5:00:00 -N cluster_${model}_${fset}_${K} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate pytorch_env > /dev/null; python ${BINDIR}/create_chmm_state_possiblelist.py main --directory ${STDIR} --outprefix $PREFIX"

