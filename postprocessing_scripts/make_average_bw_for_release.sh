#!/bin/bash
# --------------------------------------------------------
# Convert bedgraph to bw according to public data release:
# Put all of these into a directory on compbio_ce
# --------------------------------------------------------
# Run as: 
# qsub -cwd -t 1-654 -l h_vmem=30G -l h_rt=02:00:00 -N make_average_bw -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/make_average_bw_for_release.sh"
# Directories/Config:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh
start=`date +%s`
hostname -f

# Arguments:
TASK=${SGE_TASK_ID}
INFOFILE=$DBDIR/all_released_tracks_withgroup.tsv
TMP_DIR=${TMP}/bdg_bw_release_${TASK}_${RANDOM}
TMPINFO=${TMP_DIR}/tmp_info.tsv
STOREDIR=$CIDIR/groupaverages/
mkdir -p ${TMP_DIR} $STOREDIR

# Mark and group combination to run on:
MARK=$( awk -vOFS="\t" '{print $2,$4}' $INFOFILE | sort | uniq -c | sort -nr | awk -v num=${SGE_TASK_ID} 'NR==num{print $2}' )
GROUP=$( awk -vOFS="\t" '{print $2,$4}' $INFOFILE | sort | uniq -c | sort -nr | awk -v num=${SGE_TASK_ID} 'NR==num{print $3}' )
echo "[STATUS] Running for $MARK in $GROUP"

# Subset the infofile:
awk -v group="$GROUP" '$4 == group' $INFOFILE > $TMPINFO

# Make all of the individual tracks:
datasets=(imputed observed)
for dset in ${datasets[@]}; do
    echo ${dset}
    FINALBW=${STOREDIR}/average_${MARK}_${dset}_${GROUP}.bigWig
    if [[ ! -s $FINALBW ]]; then
        # ----------------------------
        # Get the chromosome averages:
        # ----------------------------
        while read chr size; do
            CPREF=${STOREDIR}/average_${MARK}_${dset}
            CFILE=${CPREF}_${GROUP}_${chr}_mean.tsv.gz
            if [[ ! -s $CFILE ]]; then
                if [[ "$dset" == "observed" ]]; then
                    cmd="python $BINDIR/calculate_mark_average_tracks.py main --mark ${MARK} --infofile ${TMPINFO} --dataset ${dset} --dir ${CONVERTED_DIR} --outprefix ${CPREF} --chrom $chr"
                else
                    cmd="python $BINDIR/calculate_mark_average_tracks.py main --mark ${MARK} --infofile ${TMPINFO} --dataset ${dset} --dir ${IMPUTED_DIR} --outprefix ${CPREF} --chrom $chr"
                fi
                echo "$cmd"
                bash -c "$cmd"
            else
                echo "[STATUS] $dset ${chr} already exists at $CFILE"
            fi
        done < $CHROMSIZES_noY

        # ---------------------
        # Make the full bigwig:
        # ---------------------
        tmpwig=${TMP_DIR}/average_${MARK}_${dset}_${GROUP}_full.tmp.wig
        rm $tmpwig
        while read chr size; do 
            echo ${chr}
            CPREF=${STOREDIR}/average_${MARK}_${dset}
            CFILE=${CPREF}_${GROUP}_${chr}_mean.tsv.gz
            # Attach this to beginning of all tracks:
            echo "" >> $tmpwig
            echo "fixedStep chrom=${chr} start=1 step=25 span=25" >> $tmpwig
            zcat $CFILE >> $tmpwig
        done < ${CHROMSIZES_noY}

        # Clip when converting because some 25 bp intervals 
        # in fixedStep may run over the tail.
        $SFTDIR/bin/wigToBigWig -clip $tmpwig ${CHROMSIZES_noY} $FINALBW
    else
        echo "[STATUS] Final file already exists:"
        ls -sh $FINALBW
    fi

done

rm -rf $TMP_DIR

end=`date +%s`
runtime=$((end-start))
echo "Finished wig to bw aggregation and conversion in $runtime seconds."
