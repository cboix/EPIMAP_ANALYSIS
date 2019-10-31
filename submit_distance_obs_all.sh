#!/bin/bash
# =================================================
# Run all GWcorr comparisons for impobs + all files
# And aggregate/process correlations
# =================================================
start=`date +%s`
hostname -f
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES

# 1. Run all of the distance calculations (chr1 only)
qsub -cwd -P compbio_lab -t 1-2051 -l h_vmem=12G -l h_rt=24:00:00 -N io_compare_all -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_obs_all.sh"


# 2. Aggregate all of the distance calculations:
cd $CIDIR/io_compare_distance
while read file; do
    echo $file
    id=${file%_all.txt}
    samp=${id%_obs_*}
    mark=${id#*_obs_}
    awk -v samp=$samp -v mark=$mark -vOFS="\t" 'BEGIN{split("",a,""); split("",ss,"");}{
    split($1,b,"_"); a[b[2]] += 1; if(a[b[2]] <= 10){ss[b[2]] += $2} 
    if(a[b[2]] == 10){print samp, mark, b[2], ss[b[2]] / a[b[2]]}}' $file > ${id}_all_proc.txt
done < <( ls *_obs_*_all.txt )

cat *_all_proc.txt | gzip -c > all_obscomp_t10.tsv.gz

# 2b. Aggregate all ranks same mark:
rm all_obscomp_ranked_compare.txt
while read file; do
    echo $file
    id=${file%_all.txt}
    samp=${id%_obs_*}
    mark=${id#*_obs_}
    awk -v samp=$samp -v mark=$mark -vOFS="\t" 'BEGIN{a = 0}$1 ~ mark && a < 25{split($1,b,"_"); a += 1; print samp, mark, b[1], $2, a}' $file >> all_obscomp_ranked_compare.txt
done < <( ls *_obs_*_all.txt )

# 3. Analyze comparisons.
R --slave -f plot_eval_sampswaps.R 


# -----------------------
# Repeat for sample diff:
# -----------------------
# 1. Run all of the distance calculations (chr1 only)
qsub -cwd -P compbio_lab -t 1-2051 -l h_vmem=10G -l h_rt=24:00:00 -N io_diff_compare_all -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_distance_diff_obs_all.sh"


# 2. Aggregate all of the distance calculations:
cd $CIDIR/iodiff_compare_distance
while read file; do
    echo $file
    id=${file%_all.txt}
    samp=${id%_diff_*}
    mark=${id#*_diff_}
    awk -v samp=$samp -v mark=$mark -vOFS="\t" 'BEGIN{split("",a,""); split("",ss,"");}{
    split($1,b,"_"); a[b[2]] += 1; if(a[b[2]] <= 10){ss[b[2]] += $2} 
    if(a[b[2]] == 10){print samp, mark, b[2], ss[b[2]] / a[b[2]]}}' $file > ${id}_all_proc.txt
done < <( ls BSS*_diff_*_all.txt )

cat *_all_proc.txt | gzip -c > all_diffcomp_t10.tsv.gz

# 2b. Aggregate all ranks same mark:
rm all_diffcomp_ranked_compare.txt
while read file; do
    echo $file
    id=${file%_all.txt}
    samp=${id%_diff_*}
    mark=${id#*_diff_}
    awk -v samp=$samp -v mark=$mark -vOFS="\t" 'BEGIN{a = 0}$1 ~ mark && a < 25{split($1,b,"_"); a += 1; print samp, mark, b[1], $2, a}' $file >> all_diffcomp_ranked_compare.txt
done < <( ls *_diff_*_all.txt )

# 3. Analyze comparisons.
R --slave -f plot_eval_diff_secondary.R 


# 4. Plot potential problematic ones:
IODIFFDIR=$CIDIR/iodiff_compare_distance
FLAGFILE=${IODIFFDIR}/iodiff_flagged.tsv
TRACKFILE=${IODIFFDIR}/iodiff_flagged_tracks.tsv

# Make markfile table:
awk '{print $1"_"$2, $3, NR}' $FLAGFILE | sort -u > $FLAGFILE.tmp
awk -vOFS="\t" '{print $1"_"$2, $3}' $SAMPLEMARK_TAB | sort -u | join - ${FLAGFILE}.tmp | sort -k4 -n | awk -vOFS="\t" '{gsub("_","\t", $1); print $1,$2, $3}' | awk -vOFS="\t" '{print $1,$2,$3,"impute_"$1"_"$2,"difference_"$1"_"$2, "impute_"$1"_"$4}'> $TRACKFILE

# Run job
qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=1:30:00 -N impdiff_flagged -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; R --slave -f ${BINDIR}/plot_impobs_tracks.R --args ${TRACKFILE} chr1 1500000 2000000"


end=`date +%s`
runtime=$((end-start))
echo "Finished obs and diff compared to all imputed sucessfully in $runtime seconds."
