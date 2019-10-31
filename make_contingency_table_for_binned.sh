#!/bin/bash
# Contigency table comparing binned datasets:
WORKDIR=$1
INFOLIST=$2
cd $WORKDIR

prefix=$( sed "${SGE_TASK_ID}q;d" ${INFOLIST} )
echo "Comparisons for ${prefix}"

CTAB=${prefix}_contingency.txt
COMPLIST=${prefix}_comparison_list.txt
awk -v pref=$prefix 'BEGIN{a=0}a==1{print $0}$0==pref{a=1}' $CHDIR/prefix_list > $COMPLIST

# Echo # of comparisons
wc $COMPLIST

# Compare: 
while read comparison
do
    echo "-- Comparing binned data of $prefix to $comparison"
    pr -mts ${prefix}_allchr.txt ${comparison}_allchr.txt > ${prefix}_cont.tmp.txt
    awk -v comp=$comparison 'BEGIN{a=0; b=0; c=0; d=0;}
    $1==0 && $2==0{a=a+1}
    $1==0 && $2==1{b=b+1}
    $1==1 && $2==0{c=c+1}
    $1==1 && $2==1{d=d+1}
    END{print comp,a,b,c,d,NR}' ${prefix}_cont.tmp.txt >> ${CTAB}
    rm ${prefix}_cont.tmp.txt
done < $COMPLIST

rm ${COMPLIST}
