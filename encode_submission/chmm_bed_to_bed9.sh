#!/bin/bash
# ----------------------------------
# ChromHMM BED to bigBed conversion:
# ----------------------------------
start=`date +%s`
hostname -f

NSTATES=18
if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $0 [-f FILE] [OPTIONS] 
    -i     File to convert [required]
    -o     Output filepath [required]
    -n     Number of states in model [one of 15, 18, 25; default: $NSTATES]"
    exit 1
fi

while getopts i:o:n: o
do      case "$o" in
    i)      export INFILE="$OPTARG";;
    o)      export OUTFILE="$OPTARG";;
    n)		export NSTATES="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

# Directories/Config:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NSTATES

# Annotation colors:
COLFILE=$DBDIR/CHMM_${MODELNAME}_colors.tsv

for INFILE in `ls BSS*segments.bed.gz`; do
    OUTFILE=${INFILE%%.bed.gz}_bed9.bed.gz
    # Conversion
    if [[ ! -s $OUTFILE ]]; then
        echo "[STATUS] Converting $INFILE to $OUTFILE with model $MODELNAME"
        zcat $INFILE | awk -vOFS="\t" -F"\t" 'FNR==NR {nam["E"$1]=$2; col["E"$1]=$3; next} {print $1, $2,$3,nam[$4], 0,".",$2,$3,col[$4]}' $COLFILE - | gzip -c > $OUTFILE
    fi
done

end=`date +%s`
runtime=$((end-start))
echo "Finished chromHMM conversion in $runtime seconds."
