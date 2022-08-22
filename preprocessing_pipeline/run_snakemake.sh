#!/bin/bash
# Grid Engine options
#$ -cwd
#$ -P compbio_lab
#$ -l h_vmem=50G 
#$ -l h_rt=48:00:00
#$ -tc 250 # Not enough space for intermediate files
#$ -j y
#$ -b y 
#$ -V 
#$ -r y 
#$ -o /broad/compbio/cboix/DEVTRAJ/db/out/ALZ
#$ -e /broad/compbio/cboix/DEVTRAJ/db/out/ALZ
# --------------------------
# Declare date and hostname:
export TODAY=$(date +%Y%m%d)
start=`date +%s`
echo "Using $( hostname -f ) on $TODAY"
# --------------------------
# For rhel7 or rhel6, locate the python libraries:
export LD_LIBRARY_PATH=/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_3.4.2/lib:/broad/software/free/Linux/redhat_7_x86_64/pkgs/python_3.5.1/lib:$LD_LIBRARY_PATH

# Inject artificial wait time
echo "[STATUS] 1. Sleeping 1-60s"
sleep $[ ( $RANDOM % 60 )  + 1 ]s
SHDIR=${SNDIR}/.snakemake/shadow

echo "[STATUS] Running command using pipeline $1: Task is $2"
ls -sha $SHDIR; rm -rf ${SHDIR}; 
snakemake --snakefile $1 -pk --rerun-incomplete --nolock $2

echo "[STATUS] 2. Sleeping 1-60s"
sleep $[ ( $RANDOM % 60 )  + 1 ]s

echo "[STATUS] Running command AGAIN using pipeline $1: Task is $2"
ls -sha $SHDIR; rm -rf ${SHDIR}; 
snakemake --snakefile $1 -pk --rerun-incomplete --nolock $2

echo "[STATUS] 3. Sleeping 1-60s"
sleep $[ ( $RANDOM % 60 )  + 1 ]s

echo "[STATUS] Running command AGAIN using pipeline $1: Task is $2"
ls -sha $SHDIR; rm -rf ${SHDIR}; 
snakemake --snakefile $1 -pk --rerun-incomplete --nolock $2

echo "Run parameters were: ID=$SGE_TASK_ID, PARAMS=\"$@\""
echo "Run successfully finished"

end=`date +%s`
runtime=$((end-start))
echo "Finished in $runtime seconds."
