#!/bin/bash
# --------------------------------------------
# Run meme-Chip enrichment for any infofile
# - Extend all regions by padding with 100bp
# - Intersect with Fasta to get sequence
# - Submit meme-chip against motif db 
# - Parse results
# - Plot results without cls (general utility)
# - Return last for plotting with modules
# --------------------------------------------
# NOTE: INFOFILE must contain absolute paths.
INFOFILE=$1
RESULTDIR=$2
mkdir -p ${RESULTDIR}

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

# Count number of files to run on:
INUM=$( wc -l ${INFOFILE} | awk '{print $1}' )
echo "[STATUS] Running MEME/FIMO submission for $INUM files out to $RESULTDIR"

# 0. Run on DB?
if [[ "0" == "1" ]]; then
    INFILE=${CORECOORD}
    DBRES=${DML_DIR}/fimo_counts_redundant
    echo "[STATUS] Submitting jobs to get FIMO counts for $INUM files in $RESULTDIR"
    qsub -cwd -t 1-358 -l h_vmem=8G -l h_rt=01:00:00 -P compbio_lab -j y -b y -V -r y -o $DBDIR/out/annotate -N getcounts_db_${UQID} "$BINDIR/meme_enrichment/run_fimo_db.sh -i $INFILE -o $DBRES -c 10000"

    if [[ ! -s $DBRES/all_fimo.counts.tsv ]] || [[ ! -s $DBRES/all_fimo_q10.counts.tsv ]]; then
        # Once this is done, aggregate all / sort 
        rm $DBRES/all_fimo.counts.tsv $DBRES/all_fimo_q10.counts.tsv
        for id in `seq 1 358`; do 
            echo "$id"
            zcat $DBRES/regions_${id}_fimo.counts.tsv.gz| awk -F"\t" -vOFS="\t" '{sub(":.*","",$1); print $1,$2}' >> $DBRES/all_fimo.counts.tsv
            zcat $DBRES/regions_${id}_fimo_q10.counts.tsv.gz| awk -F"\t" -vOFS="\t" '{sub(":.*","",$1); print $1,$2}' >> $DBRES/all_fimo_q10.counts.tsv
        done
        gzip -f $DBRES/all_fimo.counts.tsv 
        gzip -f $DBRES/all_fimo_q10.counts.tsv
        # Sort for merging
        zcat $DBRES/all_fimo.counts.tsv.gz | sort +0 -1 | gzip -c >  $DBRES/all_fimo.counts.lex.tsv.gz
        zcat $DBRES/all_fimo_q10.counts.tsv.gz | sort +0 -1 | gzip -c >  $DBRES/all_fimo_q10.counts.lex.tsv.gz
        ll -sh $DBRES/all_fimo*tsv.gz
    fi
fi

# ------------------------
# 1. MEME submission:
# ------------------------
# NOTE: ~ 4+ hrs for 70k regions - put lim at 10:00:00
echo "[STATUS] Submitting jobs to get FIMO counts for $INUM files in $RESULTDIR"
# qsub -cwd -t 1-$INUM -l h_vmem=8G -l h_rt=10:00:00 -P compbio_lab -j y -b y -V -r y -o $DBDIR/out/annotate -N getcounts_${UQID} "$BINDIR/meme_enrichment/run_fimo.sh -i ${INFOFILE} -o ${RESULTDIR}"
qsub -cwd -t 1-$INUM -l h_vmem=8G -l h_rt=01:00:00 -P compbio_lab -j y -b y -V -r y -o $DBDIR/out/annotate -N getcounts_${UQID} "$BINDIR/meme_enrichment/run_fimo.sh -i ${INFOFILE} -o ${RESULTDIR} -d 1"

# ---------------------------
# 2. Parse/plot MEME results:
# ---------------------------
# a. Aggregate raw fimo outputs:
qsub -cwd -l h_vmem=10G -l h_rt=01:00:00 -P compbio_lab -j y -b y -V -r y -o $DBDIR/out/annotate -N aggregate_fimo_${UQID} -hold_jid getcounts_${UQID} "source activate mv_env; R --slave -f $BINDIR/meme_enrichment/aggregate_raw_fimo_results.R --args ${INFOFILE} ${RESULTDIR}"

# TODO: USING DB - AGGREGATE INDPTLY TO CREATE BACKGROUND!

# b. Compute pvals in parallel (othw too slow):
echo "[STATUS] Submitting jobs to parse raw GREAT results for $INUM files in ${RESULTDIR}"
qsub -cwd -t 1-$INUM -l h_vmem=10G -l h_rt=01:00:00 -P compbio_lab -j y -b y -V -r y -o $DBDIR/out/annotate -N parse_fimo_${UQID} -hold_jid aggregate_fimo_${UQID} "source activate mv_env; R --slave -f $BINDIR/meme_enrichment/parse_fimo_results.R --args ${INFOFILE} ${RESULTDIR}"


# c. Finally, aggregate pvals and plot w/out clusters:
qsub -cwd -l h_vmem=10G -l h_rt=01:00:00 -P compbio_lab -j y -b y -V -r y -o $DBDIR/out/annotate -N aggregate_fimo_${UQID} -hold_jid parse_fimo_${UQID} "source activate mv_env; R --slave -f $BINDIR/meme_enrichment/aggregate_parsed_fimo_results.R --args ${INFOFILE} ${RESULTDIR}"


# ----------------------
# 3. Plot MEME results:
# ----------------------
# Alone, without modules:
# plotcmd="R --slave -f ${BINDIR}/plot_modules_with_great.R --args ${RESULTDIR} ${COUTPREF} '$TAGLINE' $IMGDIR"


# Return last job ID to plot with modules
echo "parse_motifs_${UQID}"
