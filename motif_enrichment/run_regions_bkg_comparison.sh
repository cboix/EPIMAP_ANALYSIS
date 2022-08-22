#!/bin/bash
# ----------------------------------------------
# Compare a set of regions with a background set
# ----------------------------------------------
start=`date +%s`
hostname -f

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

export NSHUFFLES=25
export K=8

if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [OPTIONS] 
    -i     Input regions (BEDFILE)
    -o     Output file
    -b     Input background (default BKGDHS)
    -d     Dataset to use (defaults to collated JASPAR/HOCOMOCO/Jolma2013 motifs)
    -k     KMER length (for p-value cutoff) [default: $K]
    -s     Number of shuffles to calculate [default: $NSHUFFLES]"
    exit 1
fi

while getopts i:o:b:d:k:s: o
do      case "$o" in
    i)		export INPUTFILE="$OPTARG";;
    o)		export OUTFILE="$OPTARG";;
    b)		export BKGFILE="$OPTARG";;
    d)		export ALL_PFM="$OPTARG";;
    k)		export K="$OPTARG";;
    s)      export NSHUFFLES="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

# INPUTFILE="/home/unix/cboix/data/DEVTRAJ/db/AD_scATAC/DiffPeaks100520/Ast.binomial.ADbyBraak.Diff.sort.txt"
# BKGFILE="/home/unix/cboix/data/DEVTRAJ/db/AD_scATAC/DiffPeaks100520/CelltypeSpecif.peak.Ast.fixed.txt"
# OUTFILE="/home/unix/cboix/data/DEVTRAJ/db/AD_scATAC/test.txt"

# Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
MOTIFDIR=$DBDIR/motifs
MATCHDIR=$MOTIFDIR/matches
mkdir -p $MOTIFDIR $MATCHDIR

cd $MOTIFDIR

# export TASK=${SGE_TASK_ID}
# echo "[STATUS] TASK NUMBER $TASK"
# NO TASK.
TMP_DIR=${TMP}/motif_scans_${RANDOM}
mkdir -p ${TMP_DIR}

# Genome to use:
# GENOME=$MAINDIR/genomes/hg19/hg19.fa
GENOME=${MAINDIR}/genomes/GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

# Motifs dataset:
if [[ -z "$ALL_PFM" ]]; then
    ALL_PFM=$MOTIFDIR/collated_pfm.txt
    ALL_MOT=$MOTIFDIR/collated_motifs.txt
    echo "[STATUS] Defaulting to dataset: ${ALL_PFM}"
fi

# Format the input and bkg:
if [[ "$INPUTFILE" != "" ]]; then
    awk -vOFS="\t" '$1 ~/^chr/{print $1,$2,$3,"input",NR}' $INPUTFILE > ${TMP_DIR}/input.bed
fi

if [[ "$BKGFILE" != "" ]]; then
    awk -vOFS="\t" '$1 ~/^chr/{print $1,$2,$3,"bkg", NR}' $BKGFILE > ${TMP_DIR}/bkg.bed
fi

# Merge the two for faster intersect:
BOTHFILE=${TMP_DIR}/both.bed
BOTHFA=${TMP_DIR}/both.fa
cat ${TMP_DIR}/input.bed ${TMP_DIR}/bkg.bed | sort -k1,1V -k2,2n > $BOTHFILE

# Get the genome (takes a while to extract):
bedtools getfasta -fi $GENOME -bed $BOTHFILE -fo $BOTHFA





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

MDIR=$MAINDIR/software/MOODS-python-1.9.4.1/scripts
PATH="$PATH:$MAINDIR/MOTIF_VALIDATION/bin"
SCORE=$( MotifAddPValCutoff.sh -k "$KMERS" -s "_%dmer" -m $KEEPWITHMAX $TMP_DIR/mot.txt | awk 'NR==1{print $2}' )

FINALFILE=${TMP_DIR}/matches_wshuffles_n${NSHUFFLES}_${MOTNAM}.tsv.gz

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
    cmd="python $MDIR/moods-dna.py --sep "," -s ${BOTHFA} --p-value $(awk -vK=$K 'BEGIN{print 4**(-K)}') --lo-bg 2.977e-01 2.023e-01 2.023e-01 2.977e-01 -m $PFMLIST -o ${OUTFILE}"
    echo "$cmd"
    bash -c "time $cmd"

    # Process into tsv for legibility:
    awk -F"," -vOFS="\t" '{sub(".txt","", $2); print $1, $3,$4, $2, sprintf("%0.3f",$5)}' ${OUTFILE} | gzip -c > ${OUTFILE}.gz

fi





# ---------------------------
# Count each and every motif:
# ---------------------------
# Streaming all matches through python is the most efficient way here:
MOTPREF=${MATCHDIR}/matches_wshuffles_n${NSHUFFLES}
awk '$1 ~ />/{sub(">","",$1); print $1}' $ALL_PFM > ${TMP_DIR}/motif_list.txt

while read MOTNAM; do 
    echo $MOTNAM
    FINALFILE=${MOTPREF}_${MOTNAM}.tsv.gz

    zcat $FINALFILE | awk -vOFS="\t" '{print $1,$2,$2,$4,NR}' > ${TMP_DIR}/matches.bed
    bedtools intersect -a ${TMP_DIR}/matches.bed -b ${BOTHFILE} -wa -wb | awk -vOFS="\t" '{print $4,$5,$9}' | gzip -c > ${TMP_DIR}/matched_$MOTNAM.txt.gz

    # Count up:
    zcat ${TMP_DIR}/matched_$MOTNAM.txt.gz | awk -vOFS="\t" '{print $1,$3}' | sort | uniq -c | awk -vOFS="\t" -v motif=$MOTNAM '{print $2,$3,$1, motif}' > ${TMP_DIR}/counts_${MOTNAM}.txt
    cat ${TMP_DIR}/counts_${MOTNAM}.txt >> ${OUTFILE}
done < $TMP_DIR/motif_list.txt


# Too slow:
# python $BINDIR/motif_enrichment/count_regions_bkg_motifs.py main --infile ${TMP_DIR}/input.bed --bkgfile ${TMP_DIR}/bkg.bed --outfile ${OUTFILE} --motif_list ${TMP_DIR}/motif_list.txt --motif_pref ${MOTPREF}_


rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished motif matches ($TASK) sucessfully in $runtime seconds."
