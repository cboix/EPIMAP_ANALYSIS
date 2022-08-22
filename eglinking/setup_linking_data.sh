#!/usr/bin/env bash
# Author: Benjamin T. James
# If you want data to not be re-downloaded, make DATA_DIR beforehand and populate with those files
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

# Arguments:
SRC_DIR="$BINDIR/eglinking/src"
MASTERLIST_TXT="${DMLPREF}_hg19.core.srt.txt"
MA_MATRICES="$DBDIR/markassay_matrices/"
# ENH_INDICES="$DBDIR/ENH_masterlist_indices_ind.cp.gz" # must be present somewhere
ENH_INDICES="$DBDIR/masterlist_matindices/matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_masterlist_indices.cp.gz"
DATA_DIR="$DBDIR/linking_data";
SAMPLE_DIR="$DATA_DIR/sample";

function die() { echo "$*" 1>&2; exit 1; }

[[ $- == *i* ]] || die "Not an interactive shell, rerun with bash -i"
test -d $MA_MATRICES || die "markassay_matrices dir $MA_MATRICES not found";
test -f $MASTERLIST_TXT || die "masterlist.txt must exist ($MASTERLIST_TXT)";
test -f $ENH_INDICES || die "ENH indices must exist ($ENH_INDICES)";
test -d $SRC_DIR || die "Source dir $SRC_DIR must exist and have all Python files";

# Create the appropriate conda environment:
export PATH="$SRC_DIR:$PATH";
if [[ "$( conda info --envs | grep -c "linking" )" == "0" ]]; then
    yes | conda create -n linking -c conda-forge scikit-learn numpy scipy requests wget gawk coreutils
fi
conda activate linking

mkdir -p $DATA_DIR $DATA_DIR/out
cd $DATA_DIR

function download { # sees if file exists before downloading
    basename $1 | xargs test -f && echo "$( basename $1) already exists" || wget $1;
}

if [[ ! -s $DATA_DIR/merged_log2fpkm.pc.mtx.gz ]]; then
    cp $DBDIR/public_intermediate_data_released/rnaseq_data/merged_log2fpkm.pc.mtx.gz ${DATA_DIR}
fi
if [[ ! -s $DATA_DIR/hg19.chrom.sizes ]]; then
    download 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh37_mapping/gencode.v33lift37.annotation.gtf.gz'
    download 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes'
fi

# DNase-seq non-binary matrix:
# cp ${IMPUTED_DIR}/D*/*allchr*.cp.gz $MA_MATRICES
# cp ${IMPUTED_DIR}/D*/*allchr*.cp.gz $DATA_DIR

# Turn the compressed files into MTX (submit to queue):
# TODO: Not H3K4me2
for file in `ls $MA_MATRICES/*r25_e100*merged*_csr.cp.gz`; do
    base=$( basename $file )
    attrfile=${file%_csr.*}_attr.cp.gz
    mtfile=${base%_csr.*}.mtx
    if [[ ! -f $DATA_DIR/$base ]]; then
        cp $file $DATA_DIR
        cp $attrfile $DATA_DIR
    fi
    if [[ ! -f ${mtfile}.gz ]] || [[ $file -nt ${mtfile}.gz ]]; then
        echo "Creating ${mtfile}.gz"
        qsub -cwd -P compbio_lab -l h_vmem=15G -l h_rt=0:30:00 -N process_markmat -j y -b y -V -r y -o $DATA_DIR/out "$SRC_DIR/save_mtx.py $file; gzip $mtfile"
    fi
done

# rm H3K4me2_all*.gz

# Not non-ovl:
marks=($PWD/*all_bin_on*merged.mtx.gz);
cpmarks=($PWD/*all_bin_on_*merged_csr.cp.gz);

# -------------------------------------
# Run next three processes concurrently
# Processing steps to prep annotations:
# -------------------------------------
# Assume H3K27ac has same sample list as other marks
test -f tss.txt || ${SRC_DIR}/get_tss.py "$DATA_DIR/gencode.v33lift37.annotation.gtf.gz" "merged_log2fpkm.pc.mtx.gz" > tss.txt &

mkdir -p $SAMPLE_DIR && cd $SAMPLE_DIR

awk '$0 !~ /^$/{sub(".*BSS","BSS",$0); sub("_.*","",$0); print $0}' $IMPUTED_DIR/H3K27ac/H3K27ac_all_bin_on_mixed_impobs_r25_e100_allchr_merged_names.tsv > $DATA_DIR/mark_matrix_names.txt

if [[ ! -f shared_samples.txt ]]; then
    nl "$DATA_DIR/mark_matrix_names.txt" > mark_mat_names_numbered # cannot use BASH process substitution because needs to be re-read
    zcat "$DATA_DIR/merged_log2fpkm.pc.mtx.gz" | awk '/^BSS/ { print $1 }' | sort | uniq > $SAMPLE_DIR/rnaseq_samples.txt
    cat $SAMPLE_DIR/rnaseq_samples.txt | xargs -I@ -P0 grep @ mark_mat_names_numbered | awk '{ print $2 " " $1 }' > $SAMPLE_DIR/shared_samples.txt &
fi

test -f $SAMPLE_DIR/masterlist.bed || get_masterlist.py $ENH_INDICES < $MASTERLIST_TXT > $SAMPLE_DIR/masterlist.bed &
wait;
rm mark_mat_names_numbered

# ---------------------------
# Make all enhancer data dirs
# ---------------------------
< "$DATA_DIR/mark_matrix_names.txt" xargs -n1 mkdir

mkdir -p spt
split -l 10 shared_samples.txt 'spt/x' # creates spt/xaa spt/xab ..

# Generates sample sub-dirs and rna.txt files
for file in `find spt -type f`; do
    echo "Submitting RNA and marks for chunk: $file"
    qsub -cwd -P compbio_lab -l h_vmem=15G -l h_rt=04:00:00 -N chunk_rnadata -j y -b y -V -r y -o $DATA_DIR/out "${SRC_DIR}/sample_rna_grep.py $file $DATA_DIR/merged_log2fpkm.pc.mtx.gz"
    # qsub -cwd -P compbio_lab -l h_vmem=15G -l h_rt=04:00:00 -hold_jid process_markmat -N chunk_markdata -j y -b y -V -r y -o $DATA_DIR/out "${SRC_DIR}/sample_mark_grep.py $file masterlist.bed ${cpmarks[@]}" # OLD FOR ALL MARKS
    # Run per-mark jobs:
    for markfile in `ls $DATA_DIR/*merged_csr.cp.gz`; do
        echo $markfile
        qsub -cwd -P compbio_lab -l h_vmem=15G -l h_rt=01:00:00 -hold_jid process_markmat -N chunk_markdata -j y -b y -V -r y -o $DATA_DIR/out "${SRC_DIR}/sample_mark_grep.py $file masterlist.bed $markfile"
    done
done

rm spt/x*
rmdir spt

# Generates E7.bed, E8.bed, etc files
cd $SAMPLE_DIR
< "$DATA_DIR/mark_matrix_names.txt" get_all_enhancer_regions.py E7 E8 E9 E10 E11 E15

# ----------------------------------------------------------------------
# Pre-compute the correlation of each gene with enhancers within 1.01Mb:
# ----------------------------------------------------------------------
cd $DATA_DIR

# TODO: Dense - can't load allchr necessarily?
MARK=H3K27ac
MARK=H3K4me1
MARKLIST=(H3K27ac H3K4me1 H3K4me2 H3K4me3 H3K9ac DNase-seq)
# MARKLIST=(H3K4me2 H3K4me3 H3K9ac DNase-seq)
for MARK in ${MARKLIST[@]}; do 
    echo $MARK
    MARKFILE=$DATA_DIR/${MARK}_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged_csr.cp.gz
    ZARRFILE=$DATA_DIR/${MARK}_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.zarr
    HDFILE=$DATA_DIR/${MARK}_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5
    PAIRSFILE=$DATA_DIR/pairs_ENH_ovl_df.cp.gz
    RANDFILE=$DATA_DIR/random_pairs_ENH_ovl_df.cp.gz

    # TODO: figure out how to make this work --> Not useful in this case.
    # python $BINDIR/convert_array_to_zarr.py main --infile $MARKFILE --outfile $ZARRFILE --chunk1 100000 --clevel 9
    # Make hdf5 for markfile:
    if [[ ! -s $HDFILE ]]; then
        python $BINDIR/convert_array_to_zarr.py main --infile $MARKFILE --outfile $HDFILE --make_h5py --nomake_zarr
    fi

    CLUSTER_RUN=0
    if [[ "$CLUSTER_RUN" == "1" ]]; then
        # Precompute all enhancer - gene links within 1.01 Mb: 35M corr
        qsub -cwd -P compbio_lab -l h_vmem=150G -l h_rt=12:00:00 -N precompute_$MARK -j y -b y -V -r y -o $DBDIR/out/ChromImpute "python $BINDIR/precompute_correlations.py main --exprfile $DATA_DIR/merged_log2fpkm.pc.mtx.gz --markfile $MARKFILE --tssfile $DATA_DIR/tss.txt --mlocfile $MASTERLIST_TXT --indfile $ENH_INDICES --pairsfile $PAIRSFILE --outprefix ${MARK}_precomputed --window 1010000"
        # Precompute 25 random genes per enhancer: 51M corr
        qsub -cwd -P compbio_lab -l h_vmem=150G -l h_rt=12:00:00 -N precompute_random_$MARK -j y -b y -V -r y -o $DBDIR/out/ChromImpute "python $BINDIR/precompute_correlations.py main --exprfile $DATA_DIR/merged_log2fpkm.pc.mtx.gz --markfile $MARKFILE --tssfile $DATA_DIR/tss.txt --mlocfile $MASTERLIST_TXT --indfile $ENH_INDICES --pairsfile $RANDFILE --outprefix ${MARK}_precomputed --window 1010000 --random"
    else
        CORRFILE=${MARK}_precomputed_corr.tsv.gz
        RCFILE=${MARK}_precomputed_random_corr.tsv.gz
        if [[ ! -s $CORRFILE ]]; then 
        # Precompute all enhancer - gene links within 1.01 Mb: 35M corr
            python $BINDIR/precompute_correlations.py main --exprfile $DATA_DIR/merged_log2fpkm.pc.mtx.gz --markfile $MARKFILE --tssfile $DATA_DIR/tss.txt --mlocfile $MASTERLIST_TXT --indfile $ENH_INDICES --pairsfile $PAIRSFILE --outprefix ${MARK}_precomputed --window 1010000
        fi
        if [[ ! -s $RCFILE ]]; then
            # Precompute 25 random genes per enhancer: 51M corr
            python $BINDIR/precompute_correlations.py main --exprfile $DATA_DIR/merged_log2fpkm.pc.mtx.gz --markfile $MARKFILE --tssfile $DATA_DIR/tss.txt --mlocfile $MASTERLIST_TXT --indfile $ENH_INDICES --pairsfile $RANDFILE --outprefix ${MARK}_precomputed --window 1010000 --random
        fi

        CORRHD=${MARK}_precomputed_corr.hdf5
        RCHD=${MARK}_precomputed_random_corr.hdf5
        if [[ ! -s $CORRHD ]]; then
            python $BINDIR/convert_array_to_zarr.py main --infile $CORRFILE --outfile $CORRHD --make_h5py --nomake_zarr
        fi
        if [[ ! -s $RCHD ]]; then
            python $BINDIR/convert_array_to_zarr.py main --infile $RCFILE --outfile $RCHD --make_h5py --nomake_zarr
        fi

    fi
done


# TODO: Return precomputed as hdf5:



# Check that H3K27ac correlation matches old linking/benchmarks (should be the case)
cd $DATA_DIR
MARKFILE=$DATA_DIR/H3K27ac_all_bin_on_mixed_impobs_r25_e100_allchr_merged_csr.cp.gz
PAIRSFILE=$DATA_DIR/pairs_ENH_ovl_df.cp.gz
PCDATA=$DATA_DIR/ActivePromoterEnhancerLinks_JavierrePChIC.tsv
if [[ ! -s $DATA_DIR/tss.bed ]]; then
    awk -vOFS="\t" '{print $2,$3,$3,$1,1000,$4}' $DATA_DIR/tss.txt | sort -k1V -k2n > $DATA_DIR/tss.bed
fi

if [[ ! -s $DATA_DIR/pchic_bait_loc.bed ]]; then
    zcat $PCDATA | awk -vOFS="\t" 'NR > 1{print $1,$2,$3,$4}' | sort -u | sort -k1V -k2n > ${DATA_DIR}/pchic_bait_loc.bed
fi

if [[ ! -s ${DATA_DIR}/pchic_bait_gene.tsv ]]; then
    bedtools intersect -a ${DATA_DIR}/pchic_bait_loc.bed -b ${DATA_DIR}/tss.bed -wao | awk -vOFS="\t" 'BEGIN{print "baitID", "gene"}$8 != "."{print $4, $8}' > ${DATA_DIR}/pchic_bait_gene.tsv
fi

python $BINDIR/test_calculate_correlations.py main --exprfile $DATA_DIR/merged_log2fpkm.pc.mtx.gz --markfile $MARKFILE --tssfile $DATA_DIR/tss.txt --mlocfile $MASTERLIST_TXT --indfile $ENH_INDICES --pairsfile $PAIRSFILE --outprefix H3K27ac_test


