#!/bin/bash
# =====================================================================
# Make uniform control for DNase-seq peak calling / signal generation: 
# Up-to-date as of 06/09/2017
# =====================================================================
# Directories:
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]
then
    export DBDIR=$HOME/data/EPIMAP_ANALYSIS/db
else
    export DBDIR=$HOME/EPIMAP_ANALYSIS/db
fi
export ANNDIR=$DBDIR/Annotation
export MAPDIR=$ANNDIR/umap
export DHSDIR=${DBDIR}/DNase-seq
export CHROMSIZES=$ANNDIR/hg19.chrom.sizes

# Mappability regions
export MAP36=$MAPDIR/k36.Umap.MultiTrackMappability.filtered.bed
export MAP50=$MAPDIR/k50.Umap.MultiTrackMappability.filtered.bed
export MAP100=$MAPDIR/k100.Umap.MultiTrackMappability.filtered.bed

# ==========
# Arguments:
# ==========
RLEN=$1
NREADS=$2
if [[ $# -lt 2 ]] 
then
    echo "USAGE: $(basename $0) [RLEN] [NREADS]" >&2
    echo '  [RLEN]: Read Length' >&2
    echo '  [NREADS]: Number of reads for background' >&2
    echo '  [SEED]: (OPTIONAL) Seed for shuffleBed' >&2
    echo '  [CTRLDIR]: (OPTIONAL) Output Directory' >&2
    exit 1
fi

SEED=1414 
if [[ $# -gt 2 ]]; then SEED=$3; fi

CTRLDIR=${DHSDIR}/files/WCE/tagAlign
if [[ $# -gt 3 ]]; then CTRLDIR=$4; fi
# Directories and Files:
TMPDIR=${TMP}/tmp${RANDOM}_${RANDOM}
mkdir -p $TMPDIR $CTRLDIR

# Choose mappability:
if [[ "$RLEN" == "36" ]] 
then
    MAPPABLE=$MAP36
elif [[ "$RLEN" == "50" ]]
then
    MAPPABLE=$MAP50
elif [[ "$RLEN" == "100" ]]
then
    MAPPABLE=$MAP100
else 
    echo "Don't have the read length $RLEN mappability track."
    exit 1
fi

TMPFILE=${TMPDIR}/tmp.control.tagAlign
OUTFILE=${CTRLDIR}/Uniform_BKG_CONTROL_${RLEN}_${NREADS}.tagAlign.gz

# =================================
# Script for creating control file:
# =================================
echo "Constructing $NREADS fake tags:"
bedtools random -n $NREADS -l 1 -g $CHROMSIZES -seed $SEED > $TMPFILE

# Shuffle (Note: versions 2.25 and 2.26 shuffle impossibly slowly)
echo 'Shuffling reads in mappable regions:'
shuffleBed -i $TMPFILE -incl $MAPPABLE -g $CHROMSIZES -seed $SEED | gzip -c > $OUTFILE

# Extend:
echo "Extending shuffled reads to $RLEN:"
zcat $OUTFILE | awk -v rlen=$RLEN 'BEGIN{FS="\t";OFS="\t"} $6=="+"{$3=$2+rlen; $4="N"; $5=1000; print $0} $6=="-"{$2=$3-rlen; $4="N"; $5=1000; print $0}' > $TMPFILE

# Trim:
echo 'Removing unmappable extensions:'
bedtools intersect -a ${TMPFILE} -b ${MAPPABLE} | gzip -c > ${OUTFILE}
rm -rf ${TMPFILE}

