#!/bin/bash
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh
conda activate pytorch_env 

cd $CIDIR/
mkdir -p $CIDIR/recovery_metrics

# Overall prefix + chromosome::
CPREF=$CIDIR/recovery_metrics/recovery_metrics
chr="chr19"

# Serially submit to qsub for each mark:
marks=(H3K27ac H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3 H3K36me3 H3K27me3 H3K79me2 DNase-seq ATAC-seq H4K20me1)
for mark in ${marks[@]}; do 
    echo $mark
    if [[ ! -s ${CPREF}_${mark}_${chrom}_stats.tsv.gz ]]; then
        cmd="python $BINDIR/calculate_mark_recovery_metrics.py main --mark ${mark} --infofile ${ALL_TRACKS_TAB} --dir ${IMPUTED_DIR} --outprefix ${CPREF}_${mark} --chrom $chr"
        echo "$cmd"
        qsub -cwd -l h_vmem=60G -l h_rt=2:00:00 -N recovery_${mark} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$cmd"
        # bash -c "$cmd"
    fi
done

# Collect all files:
cd $CIDIR/recovery_metrics/
tar -cvzf rc_metrics.tar.gz recovery_metrics*_stats.tsv.gz

