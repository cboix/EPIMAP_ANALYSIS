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
rm $ARGSFILE

# Enhancer and promoter states:
echo "ENH\t-e ENH -n $NUMSTATES -r 0 $MODELARGS" > $ARGSFILE
echo "PROM\t-e PROM -n $NUMSTATES -r 0 $MODELARGS" >> $ARGSFILE
# Binarized marks and assays (punctate):
echo "DNase-seq\t-e DNase-seq -n $NUMSTATES -r 0 $MODELARGS" >> $ARGSFILE
echo "H3K27ac\t-e H3K27ac -n $NUMSTATES -r 0 $MODELARGS" >> $ARGSFILE
echo "H3K4me1\t-e H3K4me1 -n $NUMSTATES -r 0 $MODELARGS" >> $ARGSFILE
echo "H3K4me2\t-e H3K4me2 -n $NUMSTATES -r 0 $MODELARGS" >> $ARGSFILE
echo "H3K4me3\t-e H3K4me2 -n $NUMSTATES -r 0 $MODELARGS" >> $ARGSFILE
echo "H3K9ac\t-e H3K9ac -n $NUMSTATES -r 0 $MODELARGS" >> $ARGSFILE
# Raw marks and assays (punctate):
# echo "DNase-seq\t-e DNase-seq -n $NUMSTATES -r 1 $MODELARGS" >> $ARGSFILE # For DNase-seq with intensities
# echo "H3K27ac\t-e H3K27ac -n $NUMSTATES -r 1 $MODELARGS" >> $ARGSFILE
# echo "H3K4me1\t-e H3K4me1 -n $NUMSTATES -r 1 $MODELARGS" >> $ARGSFILE
# echo "H3K4me2\t-e H3K4me2 -n $NUMSTATES -r 1 $MODELARGS" >> $ARGSFILE
echo "H3K4me3\t-e H3K4me3 -n $NUMSTATES -r 1 $MODELARGS" >> $ARGSFILE
# echo "H3K9ac\t-e H3K9ac -n $NUMSTATES -r 1 $MODELARGS" >> $ARGSFILE
# TODO: Need to add H3K4me2 to the mixed_impobs table to be able to run it.


RUN_NONOVL=0
RUN_DENSE=1

# Run for all of these:
while read elem SUBARGS; do
    # Overlapping index set:
    echo $SUBARGS
    JOBID=${NUMSTATES}_${elem}_${UQID}
    echo $JOBID
    # # 1. Preprocess all:
    CNUM=$( cut -f1 ${SAMPLEMARK_TAB} | sort -u | wc -l )
    qsub -cwd -t 1-$CNUM -P compbio_lab -l h_vmem=10G -l h_rt=0:30:00 -N preprocess_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t PREPROCESS $SUBARGS"

    # # # 2. Aggregate once preprocessed:
    qsub -cwd -P compbio_lab -l h_vmem=25G -l h_rt=02:00:00 -N aggregate_${JOBID} -hold_jid preprocess_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t AGGREGATE $SUBARGS"

    if [[ "$RUN_DENSE" == "1" ]]; then
        qsub -cwd -t 1-$CNUM -P compbio_lab -l h_vmem=12G -l h_rt=0:30:00 -N preprocess_dense_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t PREPROCESS $SUBARGS -z 1"
        qsub -cwd -P compbio_lab -l h_vmem=150G -l h_rt=12:00:00 -N aggregate_${JOBID} -hold_jid preprocess_dense_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t AGGREGATE $SUBARGS -z 1" # With dense
    fi

    # # ----------------------------------
    # # Also make the non-overlapping set:
    # # ----------------------------------
    if [[ "$RUN_NONOVL" == "1" ]]; then
        SUBARGS="$SUBARGS -o 1"
        JOBID=${NUMSTATES}_${elem}_nonovl_${UQID}
        echo $JOBID
        # Non-dense process + aggregate:
        CNUM=$( cut -f1 ${SAMPLEMARK_TAB} | sort -u | wc -l )
        qsub -cwd -t 1-$CNUM -P compbio_lab -l h_vmem=10G -l h_rt=0:30:00 -N preprocess_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t PREPROCESS $SUBARGS"
        qsub -cwd -P compbio_lab -l h_vmem=25G -l h_rt=02:00:00 -N aggregate_${JOBID} -hold_jid preprocess_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t AGGREGATE $SUBARGS"

        if [[ "$RUN_DENSE" == "1" ]]; then
            # Dense process + aggregate:
            qsub -cwd -t 1-$CNUM -P compbio_lab -l h_vmem=12G -l h_rt=0:30:00 -N preprocess_dense_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t PREPROCESS $SUBARGS -z 1"
        qsub -cwd -P compbio_lab -l h_vmem=50G -l h_rt=03:00:00 -N aggregate_${JOBID} -hold_jid preprocess_dense_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t AGGREGATE $SUBARGS -z 1" # With dense
        fi
    fi
done < $ARGSFILE


# For comparison of imputed and observed datasets:
if [[ "1" == "0" ]]; then
    ARGSFILE=${TMP_DIR}/args_file.tsv
    rm $ARGSFILE
    MODELARGS="-q 1 -m 2 -g $SMERGE"
    echo "H3K27ac\t-e H3K27ac -n $NUMSTATES -r 1 $MODELARGS" >> $ARGSFILE
    echo "H3K4me1\t-e H3K4me1 -n $NUMSTATES -r 1 $MODELARGS" >> $ARGSFILE
    CNUM=$( cut -f1 ${SAMPLEMARK_TAB} | sort -u | wc -l )

    # Run for all of these:
    while read elem SUBARGS; do
        echo $SUBARGS
        JOBID=${NUMSTATES}_${elem}_${UQID}
        echo $JOBID
        qsub -cwd -t 1-$CNUM -P compbio_lab -l h_vmem=10G -l h_rt=0:30:00 -N preprocess_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t PREPROCESS $SUBARGS -z 1"
        # qsub -cwd -l h_vmem=150G -l h_rt=12:00:00 -N aggregate_${JOBID} -hold_jid preprocess_dense_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t AGGREGATE $SUBARGS -z 1" # With dense
    done < $ARGSFILE

fi

# Aggregate all: 
if [[ "0" == "1" ]]; then
    cd ${DBDIR}/
    mkdir $DBDIR/markassay_matrices/ -p
    cp ${IMPUTED_DIR}/H*/*allchr*.cp.gz $DBDIR/markassay_matrices/
    cp ${IMPUTED_DIR}/D*/*allchr*.cp.gz $DBDIR/markassay_matrices/
    cp ${CHMM_FMTDIR}/DNase-seq/*allchr*.cp.gz $DBDIR/markassay_matrices/
    tar -cvzf mark_matrices.tar.gz markassay_matrices/
fi

# 3. Make matrix:
# TODO: Intersect all 3 to make matrix
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

# Create the enhancer/promoter indices:
# NOTE: Creating intersection with H3K4me1 as well as H3K27ac
if [[ "0" == "1" ]]; then
    conda activate pytorch_env 
    python ${BINDIR}/get_enh_prom_idx.py main --mixobs --nonovl --makeplots --usemark H3K27ac --savematrices --savemtx
    python ${BINDIR}/get_enh_prom_idx.py main --mixobs --nonovl --makeplots --usemark H3K4me1 --savematrices --savemtx
    python ${BINDIR}/get_enh_prom_idx.py main --mixobs --nonovl --makeplots --usemark None --savematrices --savemtx
    python ${BINDIR}/get_enh_prom_idx.py main --mixobs --makeplots --usemark H3K27ac --savematrices --savemtx
    python ${BINDIR}/get_enh_prom_idx.py main --mixobs --makeplots --usemark H3K4me1 --savematrices --savemtx
    python ${BINDIR}/get_enh_prom_idx.py main --mixobs --makeplots --usemark None --savematrices --savemtx
fi

# --------------------------------------------------------------------------------
# Compare cluster tasks (with nonovl / without) (with intersect H3K27ac and w/out)
# --------------------------------------------------------------------------------
SUBARGS="-e $ELEM -n $NUMSTATES -q 1 -m 1 -k $NCLUST -g $SMERGE -w $MWITH"
echo $JOBID
echo $SUBARGS
NUMSTATES=18
NCLUST=300
SMERGE=2
# Elements, H3K27ac int, overlaps:
ELEMLIST=(ENH PROM DNase-seq)
MWLIST=(0 1)
OVLIST=(0 1)
# conda activate pytorch_env  # Clustering
conda activate mv_env  # For gwas

# ----------------
# Non-overlapping:
# ----------------
# For clearing gwas results:
# rm -rf ${CHMM_FMTDIR}/DNase-seq/clust/*_raw/gwas_hg/
# rm -rf ${CHMM_FMTDIR}/H3K*/clust/*_raw/gwas_hg/
# rm -rf ${CALLDIR}/${MODEL}/ENH/clust/*_raw/gwas_hg/
# rm -rf ${CALLDIR}/${MODEL}/PROM/clust/*_raw/gwas_hg/

# ELEM=ENH; MWITH=1; OVLSTAT=1; # WORKS
# ELEM=ENH; MWITH=0; OVLSTAT=1; # WORKS
# ELEM=PROM; MWITH=1; OVLSTAT=1; # WORKS
# ELEM=PROM; MWITH=0; OVLSTAT=1; # WORKS
# ELEM=DNase-seq; MWITH=1; OVLSTAT=1; # DOESNT WORK MERGE - also needs ind (or just subset to any non-zero)
# ELEM=DNase-seq; MWITH=0; OVLSTAT=1; # NEEDS DNase-seq indices (or NONE)

# ELEM=ENH; MWITH=1; OVLSTAT=0; # WORKS
# ELEM=PROM; MWITH=1; OVLSTAT=0; # WORKS
# ELEM=DNase-seq; MWITH=0; OVLSTAT=0; # NEEDS DNase-seq indices (or NONE)
# ELEM=DNase-seq; MWITH=1; OVLSTAT=0; # NEEDS DNase-seq indices (or NONE)
# ELEM=DNase-seq; MWITH=0; OVLSTAT=1; # NEEDS DNase-seq indices (or NONE)

# ELEM=ENH; MWITH=0; OVLSTAT=0; # WORKS
# ELEM=PROM; MWITH=1; OVLSTAT=1; # WORKS
GWDTXT=$TMP_DIR/gwas_dirs.txt
MODEL=observed_aux_18_on_mixed_impobs_QCUT
ls ${CHMM_FMTDIR}/DNase-seq/clust/*_raw/gwas_hg/ -d > $GWDTXT
ls ${CHMM_FMTDIR}/H3K*/clust/*_raw/gwas_hg/ -d >> $GWDTXT
ls ${CALLDIR}/${MODEL}/ENH/clust/*_raw/gwas_hg/ -d >> $GWDTXT
ls ${CALLDIR}/${MODEL}/PROM/clust/*_raw/gwas_hg/ -d >> $GWDTXT

if [[ "1" == "0" ]]; then
    # Clean all of these up for re-run:
    while read dir; do
        echo $dir
        ls $dir
        rm $dir/*.Rda
        rm $dir/*.tsv.gz
        rm $dir/permuted/*.Rda
    done < $GWDTXT
fi

# Other:
# ELEMLIST=(ENH PROM DNase-seq)
# ELEMLIST=(H3K4me1 H3K9ac H3K27ac)
ELEMLIST=(ENH PROM DNase-seq H3K4me1 H3K9ac H3K27ac H3K4me3)
MWITH=0
OVLSTAT=0 # Use overlapping only
# conda activate pytorch_env
conda activate mv_env

for ELEM in ${ELEMLIST[@]}; do
    for MWITH in ${MWLIST[@]}; do 
        # for OVLSTAT in ${OVLIST}; do
            JOBID=${NUMSTATES}_${ELEM}_o${OVLSTAT}_w${MWITH}_${SMERGE}_${UQID}
            # NOTE: Run seeded (d = 1)
            echo $JOBID
            # NOTE: Only run a few jobs of each for now:
            # SUBARGS="-e $ELEM -n $NUMSTATES -q 1 -m 1 -k $NCLUST -g $SMERGE -w $MWITH -o ${OVLSTAT} -d 1"
            # qsub -cwd -P compbio_lab -t 1 -l h_vmem=20G -l h_rt=05:00:00 -N cluster_${JOBID} -hold_jid aggregate_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t CLUSTER $SUBARGS"
            # TODO: ENR on modules.

            # For extract (-a 1 = raw) the locations:
            SUBARGS="-e $ELEM -n $NUMSTATES -q 1 -m 1 -k $NCLUST -g $SMERGE -w $MWITH -o ${OVLSTAT} -a 1"
            # qsub -cwd -P compbio_lab -t 1 -l h_vmem=20G -l h_rt=03:00:00 -N extract_${JOBID} -hold_jid aggregate_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t EXTRACT $SUBARGS"
            qsub -cwd -P compbio_lab -t 1 -l h_vmem=70G -l h_rt=10:00:00 -N gwashg_${JOBID} -hold_jid extract_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t HGENR $SUBARGS"
        # done
    done
done

# -------------------------------
# Collect the outputs from HGENR:
# -------------------------------
TFILE=$TMP_DIR/hgenr_filelist.txt
rm $TFILE -f
MODEL=observed_aux_18_on_mixed_impobs_QCUT
mkdir $IMGDIR/gwas_hg -p

# Collect plots: 
while read dir; do
    echo $dir
    ls ${dir}/*stats.png >> $TFILE
    ls ${dir}/*repr.png >> $TFILE
done < $GWDTXT

while read file; do 
    echo $file
    base=${file%%/clust*}
    elem=${base##*/}
    base=${file##*/clust/cls_}
    proc=${base%%/*}
    plot=$(basename $file)
    cp $file $IMGDIR/gwas_hg/${elem}_${proc}_${plot}
done < $TFILE

cd $IMGDIR
tar -cvzf gwas_comp.tar.gz gwas_hg
echo $IMGDIR/gwas_comp.tar.gz

# Collect tables:
mkdir $DBDIR/gwas_hg_stats -p
while read dir; do
    echo $dir
    ls ${dir}/*table.tsv >> $TFILE
done < $GWDTXT

while read file; do 
    echo $file
    base=${file%%/clust*}
    elem=${base##*/}
    base=${file##*/clust/cls_}
    proc=${base%%/*}
    plot=$(basename $file)
    cp $file $DBDIR/gwas_hg_stats/${elem}_${proc}_${plot}
done < $TFILE

cd $DBDIR
tar -cvzf gwas_stats_comp.tar.gz gwas_hg_stats
echo $DBDIR/gwas_stats_comp.tar.gz



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
        # qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=05:00:00 -N cluster_${JOBID} -hold_jid aggregate_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t CLUSTER $SUBARGS"
        # TODO: Add diff clustering distances
        # TODO: RUN / DEBUG THIS
        # -----------------------------------------------
        # Run multiple runs, seeded (-d 1) minimal (-l 1)
        # -----------------------------------------------
        # SUBARGS="$SUBARGS -o 1 -d 1 -l 1"
        SUBARGS="$SUBARGS -o 1 -d 1" # -o 1 is OVERLAPS
        qsub -cwd -P compbio_lab -t 1-1000 -l h_vmem=20G -l h_rt=05:00:00 -N cluster_${JOBID} -hold_jid aggregate_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t CLUSTER $SUBARGS"
        # qsub -cwd -P compbio_lab -t 1-2 -l h_vmem=20G -l h_rt=05:00:00 -N cluster_${JOBID} -hold_jid aggregate_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t CLUSTER $SUBARGS"
        # TODO aggregate all cls results - pick top

        # ------------------------------------
        # Enrichments and downstream plotting:
        # ------------------------------------
        # qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=01:00:00 -N plotcls_${JOBID} -hold_jid cluster_${JOBID} -hold_jid aggregate_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t PLOTCLS $SUBARGS"

        # TODO: LDSC LDAK - IMPLEMENT!
        # COMMANDLIST=(GREAT MOTIF GWAS LDSC LDAK)
        # COMMANDLIST=(GREAT MOTIF GWAS FIMO) 
        # for COMMAND in ${COMMANDLIST[@]}; do 
        #     echo $COMMAND
        #     qsub -cwd -P compbio_lab -l h_vmem=5G -l h_rt=01:00:00 -N chmm_${COMMAND}_${JOBID} -hold_jid cluster_${JOBID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_chmm_to_modules.sh -t $COMMAND $SUBARGS"
        # done
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

