#!/bin/bash
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh
conda activate pytorch_env 

cd $CIDIR/coeffv

# cmd="python $BINDIR/calculate_mark_variability.py main --mark H3K27me3 --chrom chr19 --infofile $SAMPLEMARK_TAB --dir $CONVERTED_DIR/ --outprefix test_H3K27me3"
# echo "$cmd"
# bash -c "$cmd"
# cmd="python $BINDIR/calculate_mark_variability.py main --mark H3K27me3 --chrom chr19 --infofile $ALL_TRACKS_TAB --dir $IMPUTED_DIR/ --dataset imputed --outprefix imp_H3K27me3"
# echo "$cmd"
# bash -c "$cmd"


# Serialized:
chrom="chr19"
marks=(H3K27ac H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3 H3K36me3 H3K27me3 H3K79me2 DNase-seq ATAC-seq H4K20me1)
for mark in ${marks[@]}; do 
    echo $mark
    if [[ ! -s obs_${mark}_${chrom}_npeak.tsv.gz ]]; then
        cmd="python $BINDIR/calculate_mark_variability.py main --mark $mark --chrom ${chrom} --infofile $ALL_TRACKS_TAB --dir $CONVERTED_DIR/ --outprefix obs_${mark}"
        echo "$cmd"
        bash -c "$cmd"
    fi

    if [[ ! -s imp_${mark}_${chrom}_npeak.tsv.gz ]]; then
        cmd="python $BINDIR/calculate_mark_variability.py main --mark $mark --chrom ${chrom} --infofile $ALL_TRACKS_TAB --dir $IMPUTED_DIR/ --dataset imputed --outprefix imp_${mark}"
        echo "$cmd"
        bash -c "$cmd"
    fi
done

# Collect all files:
tar -cvzf obs_cv_tracks.tar.gz obs_*.tsv.gz
tar -cvzf imp_cv_tracks.tar.gz imp_*.tsv.gz

