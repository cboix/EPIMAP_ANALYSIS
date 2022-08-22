#!/bin/bash
# ---------------------------------------
# Generate a master list for a mark:
# Following similar steps to Altius/Index
# 1. Collate / Sort / Chunk the peaks 
# 2. Run peak merging for each chunk
# 3. Create a final 0/1 file as:
#    a. Sparse in MatrixMarket fmt (mtx)
#    b. Pickled CSR Matrix for python
# ---------------------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
# 0. Load in variables:
# Default: Use all imputed and observed - will separate afterwards.
INFOFILE=${ALL_TRACKS_TAB}
PREF="alltracks"
if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [OPTIONS] 
    -m     Mark [required]
    -i     Infofile [default: $INFOFILE]
    -o     Prefix [default: $PREF]"
    exit 1
fi

while getopts m:i:o: o
do      case "$o" in
    m)		MARK="$OPTARG";;
    i)		INFOFILE="$OPTARG";;
    o)		PREF="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

TMP_DIR=${TMP}/masterlist_${MARK}_${PREF}_${RANDOM}
mkdir -p ${TMP_DIR}

# -----------------------------------
# 1. Collate all corresponding files:
# -----------------------------------
AGTAB=${DBDIR}/best_worst_impobs_agreement_011119.tsv
CONCATFILE=${TMP_DIR}/concatenated_all.narrowPeak
FILELIST=${TMP_DIR}/pkfiles.tsv
rm $FILELIST $CONCATFILE
while read cell epitope fpref; do
    if [[ "$epitope" == "$MARK" ]]; then
        if [[ "$fpref" == "impute_${cell}_${epitope}" ]];then 
            PKFILE=$PKIMP_DIR/${fpref}_peaks.summit.narrowPeak.gz
        else
            # Check quality of data:
            stats=$( grep -P "${cell}\t${epitope}\t" $AGTAB | awk '{print $3}' )
            if [[ "$stats" == "LOW" ]]; then 
                echo "${fpref} is low."
                PKFILE=$PKOBS_DIR/EMPTYFILE_peaks.summit.narrowPeak.gz
            else
                PKFILE=$PKOBS_DIR/${fpref}_peaks.summit.narrowPeak.gz
            fi
        fi
        if [[ -s ${PKFILE} ]]; then
            echo "$PKFILE" >> $FILELIST
            zcat $PKFILE >> $CONCATFILE
        fi
    fi
done < ${INFOFILE}

echo "[STATUS] Concatenated $( cat $FILELIST | wc -l ) files."

# Sort concat file:
echo "[STATUS] Sorting concatenated file."
sort-bed $CONCATFILE | gzip -c > ${CONCATFILE}.gz

# Takes 6-10 min if file was appropriately sorted before.
echo "[STATUS] Chunking collated peaks into smaller files."
zcat ${CONCATFILE}.gz | awk -v minChunkSep=10000 -v minChunkSep2=4000 -v minPerChunk=10000 \
    -v maxPerChunk=100000 -v outdir=${TMP_DIR} -f ${BINDIR}/chunk_bed.awk

# Count num chunks:
CNUM=$( ls ${TMP_DIR}/chunk*.bed | wc -l )
echo "[STATUS] Will process $CNUM chunks in parallel"

qsub -cwd -t 1-$CNUM -P compbio_lab -tc 2500 -l h_vmem=8G -l h_rt=1:30:00 -N build_chunk_${MARK}_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; R --slave -f ${BINDIR}/altindex_build_list.R --args ${TMP_DIR}"

# Merge script:
qsub -cwd -P compbio_lab -l h_vmem=30G -l h_rt=5:00:00 -N merge_chunks_${MARK}_${UQID} -hold_jid build_chunk_${MARK}_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "${BINDIR}/altindex_mergepks.sh $MARK $PREF ${TMP_DIR} "

