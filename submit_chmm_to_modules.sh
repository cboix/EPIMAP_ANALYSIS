#!/bin/bash
# -----------------------------------
# Submit aggregate marks and chromhmm
# 1. Preprocess (intersect w. index)
# 2. Aggregate
# 3. Make matrix
# 4. Cluster matrix (or other)
# -----------------------------------
NUMSTATES=18
NCLUST=300
# NCLUST=400
# ELEM="DNase-seq"
# ELEM="PROM"
ELEM="ENH"
SMERGE=2
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES

# -----------------------------
# Argument sets:
TMP_DIR=${TMP}/submit_CTOM_${RANDOM}
mkdir -p ${TMP_DIR}

ARGSFILE=${TMP_DIR}/args_file.tsv
MODELARGS="-q 1 -m 1 -g $SMERGE"
echo "ENH\t-e ENH -n $NUMSTATES -r 0 $MODELARGS" > $ARGSFILE
echo "PROM\t-e PROM -n $NUMSTATES -r 0 $MODELARGS" >> $ARGSFILE
echo "DNase-seq\t-e DNase-seq -n $NUMSTATES -r 0 $MODELARGS" >> $ARGSFILE
echo "H3K27ac\t-e H3K27ac -n $NUMSTATES -r 1 $MODELARGS" >> $ARGSFILE
echo "H3K4me1\t-e H3K4me1 -n $NUMSTATES -r 1 $MODELARGS" >> $ARGSFILE
echo "H3K4me3\t-e H3K4me3 -n $NUMSTATES -r 1 $MODELARGS" >> $ARGSFILE

# Run for all of these:
while read elem SUBARGS; do
    echo $SUBARGS
    JOBID=${NUMSTATES}_${elem}_${UQID}
    echo $JOBID
    # 1. Preprocess all:
    CNUM=$( cut -f1 ${SAMPLEMARK_TAB} | sort -u | wc -l )
    qsub -cwd -t 1-$CNUM -P compbio_lab -l h_vmem=10G -l h_rt=0:30:00 -N preprocess_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t PREPROCESS $SUBARGS"

    # 2. Aggregate once preprocessed:
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=02:00:00 -N aggregate_${JOBID} -hold_jid preprocess_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t AGGREGATE $SUBARGS"
done < $ARGSFILE

# 3. Make matrix:
# A) Scale: 
# /broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromImpute/imputed/H3K27ac/
# H3K27ac_all_bin_on_mixed_impobs_r25_e100_allchr_csr.cp.gz
# H3K27ac_all_bin_on_mixed_impobs_r25_e100_allchr_attr.cp.gz
# B) Loc (given) - but could also try with 
# /broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/binarized/DNase-seq
# DNase-seq_all_bin_on_mixed_impobs_r200_e0_allchr_csr.cp.gz
# DNase-seq_all_bin_on_mixed_impobs_r200_e0_allchr_attr.cp.gz
# C) Annotation:
# /broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls/observed_aux_18_on_mixed_impobs_QCUT/ENH
# observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_r200_e0_allchr_csr.cp.gz
# observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_r200_e0_allchr_attr.cp.gz


# --------------------------------------------
# 4. Cluster matrix (more or less as follows):
# --------------------------------------------
NUMSTATES=18
NCLUST=300
ELEMLIST=(ENH PROM) # DNase-seq)
SMLIST=(2)
MWITH=1
for SMERGE in ${SMLIST[@]}; do
    for ELEM in ${ELEMLIST[@]}; do
        JOBID=${NUMSTATES}_${ELEM}_${SMERGE}_${UQID}
        SUBARGS="-e $ELEM -n $NUMSTATES -q 1 -m 1 -k $NCLUST -g $SMERGE -w $MWITH"
        echo $JOBID
        echo $SUBARGS
        qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=05:00:00 -N cluster_${JOBID} -hold_jid aggregate_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t CLUSTER $SUBARGS"

        # ------------------------------------
        # Enrichments and downstream plotting:
        # ------------------------------------
        # qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=01:00:00 -N plotcls_${JOBID} -hold_jid cluster_${JOBID} -hold_jid aggregate_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t PLOTCLS $SUBARGS"

        # COMMANDLIST=(GREAT MOTIF GWAS LDSC LDAK)
        COMMANDLIST=(GREAT MOTIF GWAS FIMO) 
        for COMMAND in ${COMMANDLIST[@]}; do 
            echo $COMMAND
            qsub -cwd -P compbio_lab -l h_vmem=5G -l h_rt=01:00:00 -N chmm_${COMMAND}_${JOBID} -hold_jid cluster_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t $COMMAND $SUBARGS"
        done
    done
done


# For running on original data only:
NUMSTATES=18
NCLUST=300
ELEMLIST=(ENH PROM)  # DNase-seq)
SMERGE=2  # Put merge at 2 for these, to collapse enh
MWITH=1
NOCLS=1
for ELEM in ${ELEMLIST[@]}; do
    JOBID=${NUMSTATES}_${ELEM}_${SMERGE}_a${NOCLS}_${UQID}
    # RUN NO CLS (-a 1 flag)
    SUBARGS="-e $ELEM -n $NUMSTATES -q 1 -m 1 -k $NCLUST -g $SMERGE -a $NOCLS -w $MWITH"
    echo $JOBID
    echo $SUBARGS
    qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=05:00:00 -N cluster_${JOBID} -hold_jid aggregate_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t EXTRACT $SUBARGS"

    # COMMANDLIST=(GREAT MOTIF GWAS LDSC LDAK)
    # COMMANDLIST=(GREAT MOTIF GWAS FIMO) 
    # for COMMAND in ${COMMANDLIST[@]}; do 
    #     echo $COMMAND
    #     qsub -cwd -P compbio_lab -l h_vmem=5G -l h_rt=01:00:00 -N chmm_${COMMAND}_${JOBID} -hold_jid cluster_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t $COMMAND $SUBARGS"
    # done
done

