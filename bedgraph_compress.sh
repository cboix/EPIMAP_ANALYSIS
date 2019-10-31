#!/bin/bash
# ===================================================
# Bedgraph sort (one chromosome at a time) + compress
# ===================================================
INPUT=$1
OUTPUT=$2
ID=$3
TMP=$4

tmpdir=${TMP}/${RANDOM}_${RANDOM}
mkdir -p $tmpdir
tmpout=$tmpdir/${ID}_full.tmp.bedgraph

# Separate by chromosome, closing handle every N steps (closing file handle is slower, but is actually consistent):
echo "Splitting $INPUT by chromosome for sorting"
awk -v tmp=$tmpdir -v id=$ID '{a=a+1; print $0 >> tmp"/"id"_"$1".bed"; if(a == 10000){close(tmp"/"id"_"$1".bed"); a=0}}' $INPUT

# Sort and append:
rm -rf $tmpout
while read file
do
    echo "Sorting $file and appending to $tmpout"
    sort -k1,1 -k2,2n $file > $file\.tmp >> $tmpout
done < <( ls ${tmpdir}/${ID}_chr*.bed | sort -k1,2V )

# Compress:
echo "Compressing to get $OUTPUT"
gzip -c $tmpout > $OUTPUT
rm -rf $tmpdir
