#!/bin/bash
# ------------------------------------------------------
# Submit plotting tracks from ALL + by metadata groups:
# For each mark (6-12) x different resolutions/areas (3)
# NOTE: PLOT IMPUTED EVEN WHEN THERE IS NO OBSERVED.
# NOTE: PLOT AVERAGE BY GROUP (OR BY DENDROGRAM CLUSTER?)
# 
# Also select random areas of genome + submit plots 
# of avg in imp/obs.
# ------------------------------------------------------
# 0. Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

# cd ${EVAL_DIR}
for task in `seq 1 33`;do
    echo $task
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:00:00 -N plot_tracks_g${task} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; R --slave -f $BINDIR/plot_impobs_groups.R --args $task chr10 0 1000000"
    # qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:00:00 -N plot_tracks_g${task} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; R --slave -f $BINDIR/plot_impobs_groups.R --args $task chr1 1000000 1250000"
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=2:00:00 -N plot_tracks_g${task} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; R --slave -f $BINDIR/plot_impobs_groups.R --args $task chr2 0 5000000"
done 

# Other Locations?
# chr19 2000000 3000000

qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=0:30:00 -N plot_tracks_g${task} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; R --slave -f $BINDIR/plot_impobs_groups.R --args $task chr10 0 1000000"

task=0
qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=3:00:00 -N plot_tracks_g${task} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; R --slave -f $BINDIR/plot_impobs_groups.R --args $task chr10 0 1000000"
qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=5:00:00 -N plot_tracks_g${task} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; R --slave -f $BINDIR/plot_impobs_groups.R --args $task chr2 0 5000000"
# qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:00:00 -N plot_tracks_g${task} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; R --slave -f $BINDIR/plot_impobs_groups.R --args $task chr1 1000000 1250000"




end=`date +%s`
runtime=$((end-start))
echo "Finished impobs submissions in $runtime seconds."
