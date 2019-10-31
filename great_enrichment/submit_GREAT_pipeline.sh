#!/bin/bash
# ------------------------
# GO Enrichment pipeline:
# 1. GREAT submission/pulldown
# 2. Parse results
# 3. Condense?
# 4. Plot results?
# ------------------------
# NOTE: INFOFILE must contain absolute paths.
INFOFILE=$1
RESULTDIR=$2
BKGFILE=$3
mkdir -p ${RESULTDIR}

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

# Count number of files to run on:
INUM=$( wc -l ${INFOFILE} | awk '{print $1}' )
echo "[STATUS] Running GREAT submission for $INUM files out to $RESULTDIR"

# ------------------------
# 1. GREAT submission:
# NOTE: One comp can only 
# submit 5 jobs at a time!
# ------------------------
echo "[STATUS] Submitting jobs to get GREAT outputs for $INUM files $RESULTDIR"
qsub -cwd -P compbio_lab -tc 2 -t 1-$INUM -l h_vmem=2G -l h_rt=00:30:00 -j y -b y -V -r y -o $DBDIR/out/annotate -N req_great_${UQID} "$BINDIR/great_enrichment/run_GREAT_submission.sh -i ${INFOFILE} -o ${RESULTDIR} -b ${BKGFILE}"
# Alternately run in series:
# for TASK in `seq 1 $INUM`; do
#     echo $TASK
#     $BINDIR/great_enrichment/run_GREAT_submission.sh -i ${INFOFILE} -o ${RESULTDIR} -t $TASK -b ${BKGFILE}
# done

# -----------------------
# 2. Parse GREAT results:
# -----------------------
echo "[STATUS] Submitting jobs to parse raw GREAT results for $INUM files in ${RESULTDIR}"
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=01:00:00 -j y -b y -V -r y -o $DBDIR/out/annotate -N parse_great_${UQID} -hold_jid req_great_${UQID} "source activate mv_env; R --slave -f $BINDIR/great_enrichment/parse_GREAT_results.R --args ${INFOFILE} ${RESULTDIR}"


# ----------------------
# 3. Plot GREAT results:
# ----------------------
# Alone, without modules:
# plotcmd="R --slave -f ${BINDIR}/plot_modules_with_great.R --args ${RESULTDIR} ${COUTPREF} '$TAGLINE' $IMGDIR"


# Return last job ID to plot with modules
echo "parse_great_${UQID}"
