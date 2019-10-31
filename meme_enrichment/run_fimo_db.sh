#!/bin/bash
# ----------
# FIMO Pipe:
# Run fimo on chunks of db - keep all counts and q10 counts.
# GOAL: 
# Get counts over 3.6M regions x 3000 motifs 
# Put into sparse matrix - aggregate according to clusterings.
# ----------
# 0. Load in variables:
TASK=${SGE_TASK_ID}
INFILE=${CORECOORD}
RESULTDIR=${DML_DIR}/fimo_counts_redundant
CHUNKSIZE=10000
if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $(basename $0) -i [INPUTFILE] -o [RESULTDIR] (OTHERARGS)
    NOTE: Set up to chunk/work with DHS Index as default
    -i  Input file [default: $INFILE] 
    -o  Output directory for results [default: $RESULTDIR]
    -c  Chunk size [default: $CHUNKSIZE]
    -t  Task (OPTIONAL)
            line from file to run. Otherwise runs SGE_TASK_ID: $SGE_TASK_ID"
    exit 1
fi

while getopts i:o:c:t: o
do      case "$o" in
    i)		INFILE="$OPTARG";;
    o)      RESULTDIR="$OPTARG";;
    c)      CHUNKSIZE="$OPTARG";;
    t)		TASK="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

echo "[STATUS] INFOFILE: $INFILE"
echo "[STATUS] RESULTDIR: $RESULTDIR"
echo "[STATUS] CHUNKSIZE: $CHUNKSIZE"
echo "[STATUS] TASK: $TASK"

# NOTE: Load after case catches variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
start=`date +%s` # For timing the pipeline

# Directories 
TMP_DIR=${TMP}/run_FIMO_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR} ${RESULTDIR}

# Files and attributes:
TMP_FILE="${TMP_DIR}/regions_${RANDOM}.bed"
TMP_FASTA="${TMP_DIR}/regions_${RANDOM}.fa"
TMP_COUNTS="${TMP_DIR}/regions_${TASK}_fimo.counts.tsv"
TMP_Q10="${TMP_DIR}/regions_${TASK}_fimo_q10.counts.tsv"
FINAL_Q10="${RESULTDIR}/regions_${TASK}_fimo_q10.counts.tsv.gz"
FINAL_ALL="${RESULTDIR}/regions_${TASK}_fimo.counts.tsv.gz"

if [[ ! -s $FINAL_ALL ]]; then
    # NOTE: Label each peak by peak NR in the original file:
    # (Otherwise, get huge results)
    if [[ $INFILE =~ \.gz$ ]]; then
        zcat ${INFILE} | awk -vOFS="\t" -v chunk=$TASK -v cs=$CHUNKSIZE 'BEGIN{start=(chunk - 1)*cs + 1; end=chunk*cs}NR >= start && NR <= end{sub("chr","",$1); print $1,$2,$3, NR}' > ${TMP_FILE}
    else
        awk -vOFS="\t" -v chunk=$TASK -v cs=$CHUNKSIZE 'BEGIN{start=(chunk - 1)*cs + 1; end=chunk*cs}NR >= start && NR <= end{sub("chr","",$1); print $1,$2,$3, NR}' $INFILE > ${TMP_FILE}
    fi

    # Turn data to fasta:
    echo "[STATUS] Getting sequences from file:"
    fastacmd="bedtools getfasta -name -fi ${REF_GENOME_FULL} -bed ${TMP_FILE} > ${TMP_FASTA}"
    echo "[STATUS] (Running) $fastacmd"
    bash -c "$fastacmd"
    ls -sh ${TMP_FASTA}

    # Locations of MEME tools/data
    MOTIFDB_DIR=$GENOMEDIR/meme_db/motif_databases/
    FIMOPATH=$SFTDIR/meme/bin/fimo
    # Motif databases to use:
    MDBLIST=(JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme 
        CIS-BP/Homo_sapiens.meme
        HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
        EUKARYOTE/jolma2013.meme)

    # Dirs for tmp files:
    mkdir -p ${RESULTDIR}/JASPAR/ $RESULTDIR/CIS-BP/ $RESULTDIR/HUMAN/ $RESULTDIR/EUKARYOTE/


    # output # of regions total
    NREGIONS=$( wc -l $TMP_FILE | awk '{print $1}' )
    echo "nregions: $NREGIONS"
    # Run fimo over all core db + aggregate:
    rm $TMP_COUNTS
    rm $TMP_Q10
    for MD in ${MDBLIST[@]}; do 
        echo "[STATUS] Intersecting with DB: $MD"
        # For all regions - collate clusters.
        MDCOUNT=$RESULTDIR/${MD}_id_${TASK}.tsv
        if  [[ ! -s $MDCOUNT ]]; then
            echo "[STATUS] Computing $MDCOUNT"
            ${FIMOPATH} --verbosity 1 -oc ${TMP_DIR} $MOTIFDB_DIR/$MD ${TMP_FASTA}
            mv ${TMP_DIR}/fimo.tsv $MDCOUNT
        fi
        # Concatenate all:
        DIRMD=$( dirname $MD ) 
        awk -F"\t" 'NR > 1 && $1 !~ /#/ && $1 != ""{print $2"_"$1,$3}' ${MDCOUNT} | sort -u | awk -vOFS="\t" -v md=$DIRMD '{print $2,$1,md}' >> ${TMP_COUNTS}
        # Aggregate counts past certain cutoff
        awk -F"\t" 'NR > 1 && $1 !~ /#/ && $1 != "" && $9 < 0.1{print $2"_"$1,$3}' ${MDCOUNT} | sort -u | awk -vOFS="\t" -v md=$DIRMD '{print $2,$1,md}' >> ${TMP_Q10}
    done

    MOTIFMAP=${RESULTDIR}/motif_id_mapping.tsv
    if [[ ! -s $MOTIFMAP ]]; then
        # Reduce names with dictionary of all motifs 
        echo "[STATUS] Making motif mapping"
        TMPMAP=$TMP_DIR/mapmd.tsv
        rm $TMPMAP
        for MD in ${MDBLIST[@]}; do 
            DIRMD=$( dirname $MD ) 
            grep "MOTIF" $MOTIFDB_DIR/$MD | awk -v md=$DIRMD '{print $3"_"$2"_"md}' >> $TMPMAP
        done
        sort -u $TMPMAP | sort +0 -1 | awk '{print $1,NR}' > ${MOTIFMAP}
    fi

    # Reduce: col1 = region_number, col2 = motif_id, sorted by region.
    # NOTE: Will be able to read as sparse matrix or intersect with join
    # NOTE: Compresses well - around 3-4x
    awk '{print $2"_"$3,$1}' ${TMP_COUNTS} | sort +0 -1 | join - ${MOTIFMAP} | awk -vOFS="\t" '{print $2, $3}' | sort -n | gzip -c > $FINAL_ALL
    awk '{print $2"_"$3,$1}' ${TMP_Q10} | sort +0 -1 | join - ${MOTIFMAP} | awk -vOFS="\t" '{print $2, $3}' | sort -n | gzip -c > $FINAL_Q10

    # Put in motif dir
    if [[ -s $FINAL_ALL ]]; then
        echo "[STATUS] Removing following intermediates:"
        ls ${RESULTDIR}/*/*meme_id_${TASK}.tsv
        rm ${RESULTDIR}/*/*meme_id_${TASK}.tsv
    fi
fi

# Clear dir:
rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished FIMO motif counting in $runtime seconds."
