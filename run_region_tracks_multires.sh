#!/bin/bash
# ------------------------------------------
# Script to make figure 1B + 1C:
# Figure 1B: Plot two random tracks 
# from each epitope of interest
# in 3 areas with different resolutions
#
# Figure 1C: Plot heatmaps of values from 
# 2000 random 25bp locations in the genome.
# ------------------------------------------
# 0. Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

export RPATH="/broad/compbio/cboix/software/miniconda2/envs/mv_env/bin/R"

# Directories:
TMP_DIR=${TMP}/region_tracks_${RANDOM}
mkdir -p ${TMP_DIR}


# Seeding from https://stackoverflow.com/a/41962458/7820599
get_seeded_random() {
    seed="$1";
    openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null;
}

seed=1

RUNALL=1
if [[ "$RUNALL" == "1" ]]; then
    eplist=( ${imp_epitopes} ATAC-seq )
    RUNPREF=track_selection_allmark
else
    eplist=( ${main_epitopes} )
    RUNPREF=track_selection
fi

# 1. Select tracks and make TSV files
SELFILE=${TMP_DIR}/${RUNPREF}.tsv
AGTAB=${DBDIR}/best_worst_impobs_agreement_011119.tsv
rm ${SELFILE}
for epitope in ${eplist[@]}; do
    # Cells with both imputed + observed tracks:
    awk -v epitope=${epitope} '$2 == epitope{print $1}' ${ALL_TRACKS_TAB} | sort | uniq -c | awk '$1 == 2{print $2}' > ${TMP_DIR}/both_io_${epitope}.tsv.tmp
    if [[ "$( grep $epitope $AGTAB | wc -l)" == "0" ]]; then
        mv ${TMP_DIR}/both_io_${epitope}.tsv.tmp ${TMP_DIR}/both_io_${epitope}.tsv
    else 
        awk -vOFS="\t" -v epitope=${epitope} '$2 == epitope && $3 != "LOW"{print $1,$3}' $AGTAB | sort | join - ${TMP_DIR}/both_io_${epitope}.tsv.tmp > ${TMP_DIR}/both_io_${epitope}.tsv
    fi
    echo "Adding 2 pairs from $epitope with $( cat ${TMP_DIR}/both_io_${epitope}.tsv | wc -l ) pairs"
    # Pick 2 tracks:
    shuf --random-source=<( get_seeded_random ${seed} ) ${TMP_DIR}/both_io_${epitope}.tsv | awk -v ep=${epitope} -vOFS="\t" 'NR <3{print $1"_"ep, ep"_"NR}' >> ${SELFILE}
done

sort -u $SELFILE > ${SELFILE}.tmp
awk -vOFS="\t" '{print $1"_"$2, $3}' $SAMPLEMARK_TAB | sort -u | join - ${SELFILE}.tmp | sort -k2 -k1 | awk -vOFS="\t" '{gsub("_","\t", $1); print $1,$2}' > $SELFILE

# 2. Select 3 regions at different resolutions:
# NOTE: SELECT REGIONS (selected mostly adhoc)

# 3. Make plots 
source activate mv_env > /dev/null; 

RUNSFILE=${TMP_DIR}/torun.tsv
# echo "chr1 1000000 1025000 TRUE" > $RUNSFILE
echo "chr19 200000 400000 FALSE" >> $RUNSFILE
echo "chr4 2000000 3500000 FALSE" >> $RUNSFILE
# Make images:
SUFFIX=${RUNPREF}_alternating.png
collate=""
IFS=' '
while read chr start end axis cap; do
    echo "Running $chr from $start to $end"
    $RPATH --slave -f ${BINDIR}/plot_impobs_tracks.R --args ${SELFILE} $chr $start $end $axis $cap 
    collate="$collate ${chr}_${start}_${end}_${SUFFIX}"
done < $RUNSFILE

# Collate images into Figure 1B
IMGDIR=${MAINDIR}/EPIMAP_ANALYSIS/img
cd $IMGDIR/impobs_tracks

NRUN=$( wc -l ${RUNSFILE} | awk '{print $1}' )
WIDTH=$( wc -l ${RUNSFILE} | awk '{print $1 * 3.25}' )
HEIGHT=$( wc -l ${SELFILE} | awk '{print $1 / 4.0}' )
# Collate with pdfnup (PDFJAM toolbox):
collate_cmd="pdfnup --nup ${NRUN}x1 $collate --outfile ${RUNPREF}_${TODAY}.pdf --papersize \"{${HEIGHT}in, ${WIDTH}in}\""

bash -c ${collate_cmd}

# PDF also in PNG (ImageMagick toolbox):
convert -density 300 ${RUNPREF}_${TODAY}.pdf -quality 100 ${RUNPREF}_${TODAY}.png

# Plot fixes:
# 3. Make clear labels of cell type

# -------------------------------------
# Figure 1C: Heatmaps of 25bp segments:
# -------------------------------------
# 1. Select all matched tracks and make TSV files
SELFILE=${TMP_DIR}/${RUNPREF}_allfiles.tsv
rm ${SELFILE}
for epitope in ${eplist[@]}; do
    # Cells with both imputed + observed tracks:
    awk -v epitope=${epitope} '$2 == epitope{print $1}' ${ALL_TRACKS_TAB} | sort | uniq -c | awk -v ep=${epitope} -vOFS="\t" '$1 == 2{print $2, ep}' >> $SELFILE
done

acclist=( CTCF SMC3 RAD21 EP300 POLR2A )
for epitope in ${acclist[@]}; do 
    cat $ACCSIMPUTATION_TAB ${EXT_SAMPMARK_TAB} | awk -v epitope=${epitope} '$2 == epitope{print $1}' | sort | uniq -c | awk -v ep=${epitope} -vOFS="\t" '$1 == 2{print $2, ep}' >> $SELFILE
done

RANDFILE=${TMP_DIR}/rand_regions.tsv
NREGIONS=5000
SNUM=$( wc -l $SELFILE | awk '{print $1}' )
OUTPREFIX=${TMP_DIR}/regions/raw_regions
mkdir -p ${TMP_DIR}/regions/

# Pick random regions:
$RPATH --slave -f ${BINDIR}/choose_rand_regions.R --args ${CHROMSIZES_noY} ${RANDFILE} ${NREGIONS}

# Extract regions for each file:
qsub -cwd -t 1-$SNUM -P compbio_lab -tc 2500 -l h_vmem=5G -l h_rt=01:00:00 -N extract_randregion_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env >> /dev/null; $RPATH --slave -f $BINDIR/extract_bw_bins.R --args ${SELFILE} ${EXT_SAMPMARK_TAB} ${RANDFILE} ${OUTPREFIX}"

# qsub -cwd -t 2055-$SNUM -P compbio_lab -tc 2500 -l h_vmem=5G -l h_rt=01:00:00 -N extract_randregion_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env >> /dev/null; $RPATH --slave -f $BINDIR/extract_bw_bins.R --args ${SELFILE} ${EXT_SAMPMARK_TAB} ${RANDFILE} ${OUTPREFIX}"


IMGPREFIX=${TMP_DIR}/rand_region_plot
# Read in all intensities, cluster, and plot all regions for figure 1c:
qsub -cwd -P compbio_lab -l h_vmem=30G -l h_rt=01:30:00 -N plot_randregion_${UQID} -hold_jid extract_randregion_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env >> /dev/null; $RPATH --slave -f $BINDIR/plot_rand_regions.R --args ${SELFILE} ${OUTPREFIX} $IMGPREFIX"


rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished figures pipeline in $runtime seconds."
