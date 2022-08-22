#!/bin/bash
# Corrects syntax errors in narrowPeak files and adds rank as column 4
# from Anshul Kundaje

if [[ $# -lt 1 ]]
    then
    echo "USAGE: $(basename $0) [narrowPeakFile]" >&2
    echo '  [narrowPeakFile]: name of narrowPeakFile' >&2
    echo '  [sortCol]: column number to use to sort (high to low)' >&2
    echo '  [addPrefix]: (OPTIONAL) Takes values Y or N. If Y then filename is prefixed to peak rank' >&2
    exit 1
fi

npkFile=$1
if [[ ! -e ${npkFile} ]]; then echo "NarrowPeak file ${npkFile} does not exist" >&2 ; exit 1; fi

sortCol=0
if [[ $# -gt 1 ]]
then
    sortCol=$2
fi

prefix=''
if [[ $# -gt 2 ]]
then
    if [[ $3 == 'Y' ]]
    then
	prefix="$(basename ${npkFile} | sed -r -e 's/^[^_]+wgEncode//g' -e 's/Rep.*$/_/g' -e 's/\..*$//g')_"
    fi
fi

if [[ ${sortCol} -eq 0 ]]
then
    if echo ${npkFile} | grep -q 'gz'
    then
	zcat ${npkFile} | awk '{printf "%s\t%d\t%d\t'"${prefix}"'Rank_%d\t%s\t%s\t%s\t%s\t%s\t%d\n",$1,$2,$3,NR,$5,$6,$7,$8,$9,$10}'
    else
	awk '{printf "%s\t%d\t%d\t'"${prefix}"'Rank_%d\t%s\t%s\t%s\t%s\t%s\t%d\n",$1,$2,$3,NR,$5,$6,$7,$8,$9,$10}' ${npkFile}
    fi
else
    if echo ${npkFile} | grep -q '\.gz'
    then
        zcat ${npkFile} | sort -k "${sortCol}nr,${sortCol}nr" | awk '{printf "%s\t%d\t%d\t'"${prefix}"'Rank_%d\t%s\t%s\t%s\t%s\t%s\t%d\n",$1,$2,$3,NR,$5,$6,$7,$8,$9,$10}'
    else
        sort -k "${sortCol}nr,${sortCol}nr" ${npkFile} | awk '{printf "%s\t%d\t%d\t'"${prefix}"'Rank_%d\t%s\t%s\t%s\t%s\t%s\t%d\n",$1,$2,$3,NR,$5,$6,$7,$8,$9,$10}'
    fi
fi
