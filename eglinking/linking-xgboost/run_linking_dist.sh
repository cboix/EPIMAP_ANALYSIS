#!/usr/bin/env bash
SRCDIR="/home/benjames/projects/linking/linking-0504/linking-with-mark/";
DATADIR="/home/benjames/projects/linking/new-model/data/";
MARKDIR="/home/benjames/data/mark/";
SAMPLE=${1};
CHROMHMM_STATE=${2};
CUTOFF=${3:-2.5};

if [ $# -le 1 ]; then
    echo "Usage: $0 sample_name chromhmm_state [cutoff=2.5]" >&2
    exit 1
fi
echo "Using ${SRCDIR}"
BASENAME=`echo "${SAMPLE}_${CHROMHMM_STATE}_${CUTOFF}" | tr ' ' '_'`;
## Usage:
##   join -j2 ~/data/Enhancer_matrix_names.txt ~/data/enh_chromhmm_states.txt | xargs -I@ "linking.sh @ 2.5"
## or:
##   parallel --eta --memfree 12G --jobs 3 './run_linking.sh {1} {2} 0' ::::  ~/data/Enhancer_matrix_names.txt :::: ~/data/enh_chromhmm_states.txt
nice $(which time) -v $SRCDIR/linking.py --sample $SAMPLE --chromhmm-state $CHROMHMM_STATE \
     -o "${BASENAME}.txt" -l "${BASENAME}.log" 2>"${BASENAME}.err" \
     --dhs-chromhmm $DATADIR/dhs_chromhmm.hdf5 \
     --pos-metadata $DATADIR/pairs_ENH_ovl_df.cp.gz \
     --neg-metadata $DATADIR/random_pairs_ENH_ovl_df.cp.gz \
     --cutoff $CUTOFF \
     --cor $DATADIR/H3K27ac_precomputed_corr.hdf5 \
           $DATADIR/H3K4me1_precomputed_corr.hdf5 \
	   $DATADIR/H3K4me2_precomputed_corr.hdf5 \
	   $DATADIR/H3K4me3_precomputed_corr.hdf5 \
	   $DATADIR/H3K9ac_precomputed_corr.hdf5 \
	   $DATADIR/DNase-seq_precomputed_corr.hdf5 \
     --rand-cor $DATADIR/H3K27ac_precomputed_random_corr.hdf5 \
                $DATADIR/H3K4me1_precomputed_random_corr.hdf5 \
		$DATADIR/H3K4me2_precomputed_random_corr.hdf5 \
		$DATADIR/H3K4me3_precomputed_random_corr.hdf5 \
		$DATADIR/H3K9ac_precomputed_random_corr.hdf5 \
		$DATADIR/DNase-seq_precomputed_random_corr.hdf5 \
     --mark ${MARKDIR}/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5 \
            ${MARKDIR}/H3K4me1_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5 \
     --distance-cutoff 100000 \

#     --furthest-mark-index-list ${MARKDIR}/H3K27ac_furthest.txt \
