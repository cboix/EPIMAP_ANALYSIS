#!/bin/bash
# --------------------------------------------
# Create list of regions for a state 
# from a STATEBYLINE directory
# outputs tsv of chmm bin locations.
# --------------------------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh


# Default is running for n25 states:
STDIR='/broad/compbio/jernst/compbio-hp/SIGNALPREDICT6/TIER2_REORDER_MODEL/states25/STATEBYLINE/'
SRDIR=${DBDIR}/state_regionlists/
PREFIX=${SRDIR}/imputed12model

python ${BINDIR}/create_chmm_state_possiblelist.py main --directory ${STDIR} --outprefix $PREFIX

# TODO: dont rerun complete chr? 
# OR: TODO: parallelize by chr (each is ~ 30m?)
qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=5:00:00 -N cluster_${model}_${fset}_${K} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate pytorch_env > /dev/null; python ${BINDIR}/create_chmm_state_possiblelist.py main --directory ${STDIR} --outprefix $PREFIX"


# -------------------------------------------------------------
# Also do this to get the possible enhancers for a FP analysis:
# -------------------------------------------------------------
STDIR="$CALLDIR/observed_aux_18_on_mixed_impobs_QCUT/STATEBYLINE/"
PREFIX=${SRDIR}/epimap_obsaux18_model

conda activate pytorch_env; 
cmd="python ${BINDIR}/create_chmm_state_possiblelist.py main --directory ${STDIR} --outprefix $PREFIX --chrom chr19 --nstate 18"
bash -c "$cmd"

# Aggregate:
SUFFIX="chr19_possible_regions.tsv"
OUTFILE=${PREFIX}_sENH_${SUFFIX}
OUTBED=${OUTFILE%%tsv}bed
rm ${OUTFILE}
numlist=( 7 8 9 10 11 15 )
for num in ${numlist[@]}; do
    echo $num
    cat ${PREFIX}_s${num}_${SUFFIX} >> $OUTFILE
done

sort -nu $OUTFILE | awk -vOFS="\t" '{print "chr19", ($1-1) * 200 + 1, $1 * 200, sprintf("chunk%07d",NR)}' > ${OUTBED}


# ------------------
# Create the matrix:
# ------------------
IDXPREF=${OUTBED%%.bed}
cp $OUTBED ${DMLPREF}_hg19.core.srt.txt
export TMP_DIR=${TMP}/CHROMHMM_chr19_AGG_${RANDOM}
mkdir -p ${TMP_DIR}

export JOBINFO=${TMP_DIR}/job_infofile.tsv

ELEMENT=H3K27ac
echo "[STATUS] Running raw signal aggregation"
export ELEMDIR=${IMPUTED_DIR}
export DATADIR=${IMPUTED_DIR}/${ELEMENT}
export FPREFIX=${ELEMENT}_all_bin_on_matched_chromhmm_enh
mkdir -p $DATADIR
# Populate temporary infofile (id, prefix, suffix)
INFOFILE=${MATCH_TRACKS_TAB}
awk -vOFS="\t" -v mark=$ELEMENT -v obsdir=$CONVERTED_DIR/ -v impdir=$IMPUTED_DIR/ '$2 == mark{if ($3 ~/impute/){print $1,"","_"$3".wig.gz", impdir} else {print $1,"","_"$3".wig.gz", obsdir}}'  $INFOFILE > $JOBINFO

# Use this bedfile to create an index:
cmd="python ${BINDIR}/aggregate_binary_txt_to_npy.py main --infofile ${JOBINFO} --dir ${ELEMDIR} --out ${FPREFIX} --noverbose --resolution 25 --extend 100 --chrom chr19 --idxpref ${IDXPREF}"

bash -c "$cmd"

# 
# ---------------
# Turn into HDF5:
# ---------------
MARKFILE=${FPREFIX}_chr19_r25_e100_csr.cp.gz
HDFILE=${FPREFIX}_chr19_r25_e100_matrix.hdf5
python $BINDIR/convert_array_to_zarr.py main --infile $MARKFILE --outfile $HDFILE --make_h5py --nomake_zarr



