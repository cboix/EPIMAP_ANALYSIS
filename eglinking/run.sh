#!/usr/bin/env bash
# Author: Benjamin T. James
# If you want data to not be re-downloaded, make DATA_DIR beforehand and populate with those files
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

# Arguments:
SRC_DIR="$BINDIR/eglinking/src"
MASTERLIST_TXT="${DMLPREF}_hg19.core.srt.txt"
MA_MATRICES="$DBDIR/markassay_matrices/"
ENH_INDICES="$DBDIR/ENH_masterlist_indices_ind.cp.gz" # must be present somewhere

# Directories for linking:
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
yes | conda create -n linking -c conda-forge scikit-learn numpy scipy requests wget gawk coreutils
conda activate linking

mkdir -p $DATA_DIR $DATA_DIR/out
cd $DATA_DIR

function download { # sees if file exists before downloading
    basename $1 | xargs test -f && echo "$( basename $1) already exists" || wget $1;
}

cp $DBDIR/public_intermediate_data_released/rnaseq_data/merged_log2fpkm.pc.mtx.gz ${DATA_DIR}
cp $DBDIR/public_intermediate_data_released/enhancers_data/Enhancer_matrix_names.txt ${DATA_DIR}
download 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh37_mapping/gencode.v33lift37.annotation.gtf.gz'
download 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes'

# Turn the compressed files into MTX (submit to queue):
for file in `ls $MA_MATRICES/*merged*_csr.cp.gz`; do
    base=$( basename $file )
    mtfile=${base%_csr.*}.mtx
    if [[ ! -f ${mtfile}.gz ]] || [[ $file -nt ${mtfile}.gz ]]; then
        echo "Creating ${mtfile}.gz"
        qsub -cwd -P compbio_lab -l h_vmem=15G -l h_rt=0:30:00 -N process_markmat -j y -b y -V -r y -o $DATA_DIR/out "$SRC_DIR/save_mtx.py $file; gzip $mtfile"
    fi
done

marks=($PWD/*merged.mtx.gz);

# -------------------------------------
# Run next three processes concurrently
# Processing steps to prep annotations:
# -------------------------------------
# Assume H3K27ac has same sample list as other marks
test -f tss.txt || ${SRC_DIR}/get_tss.py "$DATA_DIR/gencode.v33lift37.annotation.gtf.gz" "merged_log2fpkm.pc.mtx.gz" > tss.txt &

mkdir -p $SAMPLE_DIR && cd $SAMPLE_DIR

nl "$DATA_DIR/Enhancer_matrix_names.txt" > enh_mat_names_numbered # cannot use BASH process substitution because needs to be re-read
zcat "$DATA_DIR/merged_log2fpkm.pc.mtx.gz" | awk '/^BSS/ { print $1 }' | sort | uniq > rnaseq_samples.txt
cat rnaseq_samples.txt | xargs -I@ -P0 grep @ enh_mat_names_numbered | awk '{ print $2 " " $1 }' > shared_samples.txt &

test -f masterlist.bed || get_masterlist.py $ENH_INDICES < $MASTERLIST_TXT > masterlist.bed &
wait;
rm enh_mat_names_numbered

# ---------------------------
# Make all enhancer data dirs
# ---------------------------
< "$DATA_DIR/Enhancer_matrix_names.txt" xargs -n1 mkdir

mkdir -p spt
split -l 10 shared_samples.txt 'spt/x' # creates spt/xaa spt/xab ..

# Generates sample sub-dirs and rna.txt files
for file in `find spt -type f`; do
    echo "Chunk: $file"

    sample_rna_grep.py

    find spt -type f | xargs -P0 -I@ sample_rna_grep.py @ "$DATA_DIR/merged_log2fpkm.pc.mtx.gz"

    # find spt -type f | xargs -P0 -I@ sample_mark_grep.py @ masterlist.bed ${marks[@]} # Will take a while

done

find spt -type f | xargs -P0 -I@ sample_rna_grep.py @ "$DATA_DIR/merged_log2fpkm.pc.mtx.gz"
find spt -type f | xargs -P0 -I@ sample_mark_grep.py @ masterlist.bed ${marks[@]} # Will take a while


rm spt/x*
rmdir spt


# Generates E7.bed, E8.bed, etc files
< "$DATA_DIR/Enhancer_matrix_names.txt" get_all_enhancer_regions.py E7 E8 E9 E10 E11 E15
