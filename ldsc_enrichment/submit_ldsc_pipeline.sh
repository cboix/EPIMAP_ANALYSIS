#!/bin/bash
# ------------------------------------
# LDSC Functional Enrichment pipeline:
# 1. Turn data to annot and l2 files
# 2. Run regressions for each GWAS
# 3. Parse + plot results for LDSC
# ------------------------------------
# NOTES:
# Must have pre-processed GWAS into sumstats
# Var INFOFILE must contain absolute paths.
# For now, all INFOFILE files should all 
#   be located in the same folder
# ------------------------------------
INFOFILE=$1
RESULTDIR=$2
mkdir -p ${RESULTDIR}

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

# ---------------------------------------------------
# 1. Process annotations: generate annot and l2 files
# ---------------------------------------------------
# Count number of files to run on:
INUM=$( wc -l ${INFOFILE} | awk '{print $1}' )
echo "[STATUS] Running LDSC processing for $INUM files"

# Submit processing job:
qsub -cwd -t 1-$INUM -P compbio_lab -l h_vmem=5G -l h_rt=1:00:00 -N proc_bedannot_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "${BINDIR}/ldsc_enrichment/process_bed_to_annot.sh $INFOFILE"

# ------------------------------------------------------------
# 2. Partition heritability command for each GWAS of interest:
# ------------------------------------------------------------
# Count number of GWAS:
export SSLIST=${DBDIR}/sumstats_gwas_list.tsv
export MAINGWAS=/broad/compbio/data/gwas
ls $MAINGWAS/*/sumstats/*.sumstats.gz | awk -vOFS="\t" '{a=$1; split($1,arr,"/"); sub(".sumstats.gz","",arr[8]); print $1, arr[6], arr[8]}' > $SSLIST
GNUM=$( cat $SSLIST | wc -l )

# Submit partition job:
qsub -cwd -t 1-$GNUM -P compbio_lab -l h_vmem=10G -l h_rt=1:00:00 -N part_h2_${UQID} -hold_jid proc_bedannot_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "${BINDIR}/ldsc_enrichment/run_partition_heritability_per_gwas.sh $INFOFILE $SSLIST $RESULTDIR"

# for above:
# python ldsc.py 
# 	--h2 BMI.sumstats.gz\
# 	--ref-ld-chr baseline.\ 
# 	--w-ld-chr weights.\
# 	--overlap-annot\
# 	--frqfile-chr 1000G.mac5eur.\
# 	--out BMI_baseline

# ------------------------------------------------------
# 3. Parse and plot results of partitioning heritability
# ------------------------------------------------------
# TODO: Decide on opts etc.
qsub -cwd -t 1-$INUM -P compbio_lab -l h_vmem=10G -l h_rt=1:00:00 -N proc_gwas_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "${BINDIR}/ldsc_enrichment/plot_LDSC_results.sh $INFOFILE $GWASLIST $RESULTDIR"


# ------------------------------------------
# 4. Alternative run scheme (cell-specific):
# ------------------------------------------
# python ldsc.py 
# 	--h2 BMI.sumstats.gz\
# 	--w-ld-chr weights.\
# 	--ref-ld-chr CNS.,baseline.\
# 	--overlap-annot\
# 	--frqfile-chr 1000G.mac5eur.\
# 	--out BMI_CNS\
# 	--print-coefficients

