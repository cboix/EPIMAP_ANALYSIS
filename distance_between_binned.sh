#!/bin/bash
# ====================================================
# Compares sets of CHMM binned datasets to each other.
# Returns:
# 1. Jaccard distance
# 2. Mutual information
# 3. Core sets for the same epitope
# 4. Conditional mutual information for same epitope.
# ====================================================
export PROJECT=$1
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]
then
    export DBDIR=$HOME/data/${PROJECT}/db
    export BINDIR=$HOME/data/EPIMAP_ANALYSIS/bin
else
    export DBDIR=$HOME/${PROJECT}/db
    export BINDIR=$HOME/EPIMAP_ANALYSIS/bin
fi
export CHDIR=$DBDIR/ChromHMM/binarized

# Catalog available data:
# -----------------------
cd $CHDIR
ls -d */ > $CHDIR/epitope_list
ls -d */*_chr1_binary.txt | awk '{sub("_chr1_binary.txt","",$0); print $0}' | sort -u > $CHDIR/prefix_list

# Collate chrs:
# -------------
IFS=$'\t'
while read prefix
do 
    echo $prefix
    # NOTE: numerical ordering does not matter (don't compare on Y)
    cat ${prefix}_chr[0-9X]*_binary.txt > ${prefix}_allchr.txt
done < $CHDIR/prefix_list

# Make contingency tables:
# ------------------------
# Submit as array job:
PREFNUM=$( awk 'END{print NR}' $CHDIR/prefix_list )
qsub -cwd -t 1-$PREFNUM -l h_vmem=20G -N ${PROJECT}_contingency -j y -b y -V -r y -o $DBDIR/out "$BINDIR/make_contingency_table_for_binned.sh $CHDIR $CHDIR/prefix_list" 

CONT=$CHDIR/full_contingency_table.txt
rm -rf $CONT
while read prefix
do
    awk -v pref=$prefix '{print pref,$0}' ${prefix}_contigency.txt >> $CONT
done < $CHDIR/prefix_list

# Once done, calculate distances:
# -------------------------------
R --slave -f $BINDIR/distance_from_contingency.R $CONT
 
# 1. Jaccard distance:
# -- Intersection/Union

# 2. Mutual information:
# -- sum p(x,y) log(p(x,y)/(p(x)p(y))) for (x,y) in (0,0), (0,1), (1,0), (1,1).


# Comparisons within each of the epitope types:
# for 5? 
while read epitope 
do 
    echo "$epitope"
    grep $epitope $CHDIR/prefix_list | awk '{sub("/","\t",$0); print $2}' | sort | > ${epitope}cell_types

    FULL=${epitope}full_allchr.txt
    rm -rf $FULL
    while read ct
    do
        if [[ ! -s $FULL ]] 
        then
            cp ${epitope}${ct}_allchr.txt ${FULL}
        else
            pr -mts $FULL ${epitope}${ct}_allchr.txt > ${FULL}.tmp
            mv ${FULL}.tmp $FULL
        fi
    done < ${epitope}cell_types
    # Remove awkward headers:
    awk '$1 == 1 || $1 == 0' $FULL > ${FULL}.tmp
    mv ${FULL}.tmp $FULL

    Rscript -e '
    #
    # 3. Core sets for the same epitope
    # 4. Conditional mutual information for same epitope.

done < $CHDIR/epitope_list


