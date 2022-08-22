#!/bin/bash
# get and process DNase data:
export DNASE=$DBDIR/DNase
DNFILES=$ANNDIR/DNase_files.txt

if [[ !-s ${DNASE}/ENCFF001SOG\.bed\.gz ]]
then
    cd $DNASE
    xargs -n 1 curl -O -L < $( awk 'NR==1{print $0}$0 ~ /bed\.gz/{print $0}' $DNFILES )
fi

# DNase Metadata:
awk 'NR==1{print $0}$2 == "bed"{print $0}' ${DNASE}/metadata.tsv > ${DNASE}/bed_metadata.tsv
# Determine useful data:
# NOTE some is broadpeak, some narrowpeak

# Normalize individual files:


# Merge Replicates


