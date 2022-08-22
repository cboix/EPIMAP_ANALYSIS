#!/bin/bash
#SBATCH -J score_IC_tracks        # Job name (-N)
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G              # Job memory request (-l h_vmem)
#SBATCH --time=01:00:00         # Time limit hrs:min:sec (-l h_rt)
#SBATCH -p kellis
#SBATCH --output=/home/cboix/data/EPIMAP_ANALYSIS/db/out/IC/score_IC_tracks_%j.log
#SBATCH --export=ALL            # Copy env. variables (-V)
#SBATCH --requeue               # For if requeue needed (-r y)

# Grid Engine options
#$ -N score_IC_tracks
#$ -cwd
# #$ -P compbio_lab
#$ -l h_vmem=10G 
#$ -l h_rt=00:20:00
#$ -tc 250
#$ -j y
#$ -b y 
#$ -V 
#$ -r y 
#$ -o /broad/compbio/cboix/EPIMAP_ANALYSIS/db/out/IC
#$ -e /broad/compbio/cboix/EPIMAP_ANALYSIS/db/out/IC
#$ -t 1-1326

# 1326 runs:
# sbatch --array=1-1326 $BINDIR/ic_eval/run_IC_eval.sh
# qsub -t 1-1326 -hold_jid process_IC_files $BINDIR/ic_eval/run_IC_eval.sh 

# Config and directories:
# -----------------------
source ${HOME}/data/EPIMAP_ANALYSIS/bin/ic_eval/config_IC_files.sh

# Cluster type: 
if [[ "$SLURM_JOB_ID" == "" ]]; then
    TASK=$SGE_TASK_ID 
else 
    TASK=$SLURM_ARRAY_TASK_ID
fi

# ------------------------------------------------
# Figure out which track to process for each task:
# ------------------------------------------------
nteam=$( wc -l $TEAMFILE | cut -f1 -d" " )
ntrack=$( wc -l $INFOFILE | cut -f1 -d" " )
let "ngrp=$nteam"

# Track:
let "tracknum=$TASK % $ntrack"
sample=$( awk -v id=$tracknum 'NR==(id+1){print $1}' ${INFOFILE} )
mark=$( awk -v id=$tracknum 'NR==(id+1){print $2}' ${INFOFILE} )
track=$( awk -v id=$tracknum 'NR==(id+1){print $3}' ${INFOFILE} )

# Group:
let "grpnum=($TASK - $tracknum) / $ntrack + 1"  # Teams 1-26
num=$( awk -v id=$grpnum 'NR==id{print $1}' ${TEAMFILE} )
team=$( awk -v id=$grpnum 'NR==id{print $2}' ${TEAMFILE} )

cd $ICDIR

# Compare the imputed files with the blind data:
# ----------------------------------------------
export EVALFILE=${STATDIR}/${sample}${mark}_${num}_eval.tsv
export BLINDBASE=blind_${sample}${mark}.bedgraph.wig.gz
export BASEWIG=${num}_${track%%.*}.bedgraph.wig.gz
echo "$team -- $num -- $sample $mark"

if [[ ! -s $EVALFILE ]]; then
    echo "[STATUS] Running scoring for $team (# $num): track $BASEWIG"
    cmd="$CIBASE Eval $WIGDIR $BLINDBASE $WIGDIR $BASEWIG $CHROMFILE > ${EVALFILE}"
    echo "$cmd"
    bash -c "$cmd"
fi

# Evaluation for qnorm data:
# --------------------------
export QEVALFILE=${STATDIR}/${sample}${mark}_${num}_evalqnorm.tsv
export QBLINDBASE=qnorm_${sample}${mark}.bedgraph.wig.gz
if [[ ! -s $QEVALFILE ]]; then
    echo "[STATUS] Running scoring for $team (# $num): track $BASEWIG"
    cmd="$CIBASE Eval $WIGDIR $QBLINDBASE $WIGDIR $BASEWIG $CHROMFILE > ${QEVALFILE}"
    echo "$cmd"
    bash -c "$cmd"
fi

# TODO: Evaluate with the two blacklists

end=`date +%s`
runtime=$((end-start))
echo "Finished track evaluation sucessfully in $runtime seconds."
