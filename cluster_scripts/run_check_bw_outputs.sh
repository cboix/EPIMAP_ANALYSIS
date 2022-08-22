# Grid Engine options
#$ -N check_bwout
#$ -cwd
#$ -P compbio_lab
#$ -l h_vmem=10G 
#$ -l h_rt=01:00:00  # Some of the BQ datasets are very large.
#$ -tc 250
# #$ -M cboix@mit.edu 
# #$ -m a 
#$ -j y
#$ -b y 
#$ -V 
#$ -r y 
#$ -o /broad/compbio/cboix/EPIMAP_ANALYSIS/db/out/ChromImpute
#$ -e /broad/compbio/cboix/EPIMAP_ANALYSIS/db/out/ChromImpute
#$ -t 1-869

# -------------------------------------------------------------
# Check all converted files and imputed outputs for each sample
# Run as: qsub "$BINDIR/run_check_bw_outputs.sh"
# -------------------------------------------------------------
start=`date +%s`
hostname -f

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES

SAMPLE=$( cut -f1 ${SAMPLEMARK_TAB} | sort -u | sed "${SGE_TASK_ID}q;d" )

echo "[STATUS] Running on $SAMPLE"

# ---------------
# Converted files
# ---------------
for cfile in `ls $CONVERTED_DIR/chr[0-9X]*_FINAL_*_${SAMPLE}.sub*wig.gz`; do 
    NL=$(zcat $cfile | uniq | wc -l )
    if [[ "$NL" == "3" ]]; then
        bn=$( basename $cfile )
        chr=${bn%%_FINAL*}
        chrlen=${#chr}
        # echo $chrlen
        if [ $chrlen -lt 6 ]; then
            base=${cfile#*_FINAL_}
            mark=${base%%_*}
            echo "[FLAGGED] $SAMPLE $mark $chr has $NL uniq"
            echo "$SAMPLE\t$mark\t$chr\t$cfile" >> $CIDIR/flag_converted_tracks.txt
        fi
    fi
done


# Imputed files:
for cfile in `ls $IMPUTED_DIR/chr[0-9X]*_impute_${SAMPLE}_*wig.gz`; do 
    NL=$(zcat $cfile | uniq | wc -l )
    if [[ "$NL" == "3" ]]; then
        bn=$( basename $cfile )
        chr=${bn%%_FINAL*}
        chrlen=${#chr}
        # echo $chrlen
        if [ $chrlen -lt 6 ]; then
            base=${cfile#*_FINAL_}
            mark=${base%%_*}
            echo "[FLAGGED] $SAMPLE $mark $chr has $NL uniq"
            echo "$SAMPLE\t$mark\t$chr\t$cfile" >> $CIDIR/flag_imputed_tracks.txt
        fi
    fi
done


# Check results for converted tracks:
echo "[STATUS] Tracks for converted:"
grep "^$SAMPLE" $CIDIR/flag_converted_tracks.txt | awk -vOFS="\t" '{print $1,$2,$3}' | sort -k3V

while read samp mark; do
    echo "[STATUS] Looking for original $mark data"
    if [[ "$mark" == "DNase-seq" ]] || [[ "$mark" == "ATAC-seq" ]]; then
        TADIR=$DBDIR/$mark/files/$mark/tagAlign
    else 
        TADIR=$DBDIR/ChIP-seq/files/$mark/tagAlign
    fi
    # Look at orig track  loc:
    ls -sh $TADIR/*$SAMPLE*
    tafile=$TADIR/FINAL_${mark}_${SAMPLE}.sub.tagAlign.gz 
    if [[ -s $tafile ]]; then
        echo "[STATUS] ${tafile} has"
        zcat $TADIR/FINAL_${mark}_${SAMPLE}.sub.tagAlign.gz | cut -f1 | uniq -c 
    fi
done < <( grep "^$SAMPLE" $CIDIR/flag_converted_tracks.txt | awk -vOFS="\t" '{print $1,$2}' | sort -u )

# For imputed tracks
echo "[STATUS] Tracks for imputed:"
grep "^$SAMPLE" $CIDIR/flag_imputed_tracks.txt | awk -vOFS="\t" '{print $1,$2,$3}' | sort -k3V

# Report number of issues:
NCONV=$( grep "^$SAMPLE" $CIDIR/flag_converted_tracks.txt | wc -l) 
NIMP=$( grep "^$SAMPLE" $CIDIR/flag_imputed_tracks.txt | wc -l )
echo "[STATUS] Found $NCONV converted issues and $NIMP imputed issues"

end=`date +%s`
runtime=$((end-start))
echo "[STATUS] Finished check sample $SAMPLE sucessfully in $runtime seconds."
