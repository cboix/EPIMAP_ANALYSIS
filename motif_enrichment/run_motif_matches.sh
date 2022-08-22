#!/bin/bash
# -------------------------------------------------------
# Compute motif (and shuffles) matches across the genome:
# -------------------------------------------------------
start=`date +%s`
hostname -f

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

export NSHUFFLES=100
export NKEEP=10
export K=8
export GEN="hg19"

if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [OPTIONS] 
    -d     Dataset to use (defaults to collated JASPAR/HOCOMOCO/Jolma2013 motifs)
    -k     KMER length (for p-value cutoff) [default: $K]
    -g     Genome to use (hg19/hg38) [default: $GEN]
    -s     Number of shuffles to calculate [default: $NSHUFFLES]"
    exit 1
fi

while getopts d:k:g:s: o
do      case "$o" in
    d)		export ALL_PFM="$OPTARG";;
    k)		export K="$OPTARG";;
    g)		export GEN="$OPTARG";;
    s)      export NSHUFFLES="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done


# Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
MOTIFDIR=$DBDIR/motifs
MATCHDIR=$MOTIFDIR/matches
mkdir -p $MOTIFDIR $MATCHDIR

export TASK=${SGE_TASK_ID}
echo "[STATUS] TASK NUMBER $TASK"

# Genome to use:
if [[ "$GEN" == "hg19" ]]; then
    GENOME=$MAINDIR/genomes/hg19/hg19.fa
elif [[ "$GEN" == "hg38" ]]; then
    GENOME=${MAINDIR}/genomes/GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
else 
    echo "[EXITING] Genome ($GEN) not recognized, pick either hg19 or hg38"
    exit 1
fi

# Motifs dataset:
if [[ -z "$ALL_PFM" ]]; then
    ALL_PFM=$MOTIFDIR/collated_pfm.txt
    ALL_MOT=$MOTIFDIR/collated_motifs.txt
    echo "[STATUS] Defaulting to dataset: ${ALL_PFM}"
fi

# Scripts for MOODS motif scanning:
MDIR=$MAINDIR/software/MOODS-python-1.9.4.1/scripts
PATH="$PATH:$MAINDIR/MOTIF_VALIDATION/bin"  # For some utils; don't need anymore

# Actual arguments to the datasets:
TASK=${SGE_TASK_ID}
TMP_DIR=${TMP}/motif_scans_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR}


# --------------
# Extract motif:
# --------------
awk -v num=$TASK '$1 ~ />/{a=a+1}{if(a == num){print $0}}' $ALL_PFM > ${TMP_DIR}/pfm.txt
MOTNAM=$( awk '$1 ~ />/{sub(">","",$1); print $1}' ${TMP_DIR}/pfm.txt )
echo "Motif is: $MOTNAM"
echo "PFM is:"
cat ${TMP_DIR}/pfm.txt

# Get motif:
awk -v mot=$MOTNAM '$1 ~ />/{a = 0}$1 == ">"mot{a=1}{if(a == 1){print $0}}' $ALL_MOT > ${TMP_DIR}/mot.txt

KMERS=8
KEEPWITHMAX=0
PATH="$PATH:$MAINDIR/MOTIF_VALIDATION/bin"
SCORE=$( MotifAddPValCutoff.sh -k "$KMERS" -s "_%dmer" -m $KEEPWITHMAX $TMP_DIR/mot.txt | awk 'NR==1{print $2}' )

FINALFILE=${MATCHDIR}/matches_wshuffles_n${NSHUFFLES}_${MOTNAM}.tsv.gz
if [[ "$GEN" == "hg38" ]]; then
    FINALFILE=${MATCHDIR}/matches_${GEN}_wshuffles_n${NSHUFFLES}_${MOTNAM}.tsv.gz
fi

if [[ ! -s $FINALFILE ]]; then
    # ---------------------------
    # Create the shuffled motifs:
    # ---------------------------
    echo "[STATUS] Making shuffles"
    python $BINDIR/motif_enrichment/motif_shuffle.py main --infile $TMP_DIR/pfm.txt --outprefix ${TMP_DIR}/pfm --nshuffle $NSHUFFLES

    # Create pfm list:
    PFMLIST="$TMP_DIR/pfm.txt"
    for file in `ls $TMP_DIR/pfm_*.txt`; do
        PFMLIST="$PFMLIST $file"
    done

    # --------------------------------------
    # Run full set of matches on the genome:
    # --------------------------------------
    # Note: 1 + 25 takes 4m; 1 + 100 takes ~14 minutes
    echo "[STATUS] Running matches"
    OUTFILE=$TMP_DIR/scan_wshuffles_n${NSHUFFLES}.tsv
    cmd="python $MDIR/moods-dna.py --sep "," -s ${GENOME} --p-value $(awk -vK=$K 'BEGIN{print 4**(-K)}') --lo-bg 2.977e-01 2.023e-01 2.023e-01 2.977e-01 -m $PFMLIST -o ${OUTFILE}"
    echo "$cmd"
    bash -c "time $cmd"

    # Process into tsv for legibility:
    awk -F"," -vOFS="\t" '{sub(".txt","", $2); print $1, $3,$4, $2, sprintf("%0.3f",$5)}' ${OUTFILE} | gzip -c > ${OUTFILE}.gz

    # Move to matchdir:
    mv ${OUTFILE}.gz $FINALFILE

    # TODO: Keep specific shuffles?
else
    echo "[STATUS] File already exists at: $FINALFILE"
    ls -sh $FINALFILE
fi


rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished motif matches ($TASK) sucessfully in $runtime seconds."
