#!/bin/bash
# Aggregate eval statistics for IC data:

source ${HOME}/data/EPIMAP_ANALYSIS/bin/ic_eval/config_IC_files.sh

cd $STATDIR

FULLFILE=$STATDIR/aggregated_IC_stats_simple.tsv
QFULLFILE=$STATDIR/aggregated_IC_stats_qnorm.tsv
awk -vOFS="\t" 'NR==1{print $0,"team","track"}' $STATDIR/C05M17_100_eval.tsv > $FULLFILE
awk -vOFS="\t" 'NR==1{print $0,"team","track"}' $STATDIR/C05M17_100_evalqnorm.tsv > $QFULLFILE
while read file; do 
    base=${file%_*}
    track=${base%_*}
    team=${base#*_}
    echo $team $track
    QEVAL=$STATDIR/${base}_evalqnorm.tsv
    awk -vOFS="\t" -v team=$team -v track=$track 'NR==2{print $0,team, track}' $file | awk '$0 !~ /Native/' >> $FULLFILE
    awk -vOFS="\t" -v team=$team -v track=$track 'NR==2{print $0,team, track}' $QEVAL | awk '$0 !~ /Native/' >> $QFULLFILE
done < <(ls *_eval.tsv)


rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished track download and conversion sucessfully in $runtime seconds."
