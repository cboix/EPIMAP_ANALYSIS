#!/bin/bash
# Put together all ChromImpute distance measures:
cd ${DISTANCE_DIR}
# cd $IMPDIST_DIR

rm ${DISTANCE_DIR}/distance_all_*.tsv

# MERGE:
while read mark
do
    echo $mark
    FULL=${DISTANCE_DIR}/distance_all_${mark}.tsv

    # List available:
    ls ${DISTANCE_DIR}/BSS*_${mark}.txt | awk -v mark=$mark '{sub(".*/","",$0); sub("_[A-Za-z0-9-]*.txt","",$0); print $0}' | sort +0 -1 > $CIDIR/template_${mark}

    # REMEMBER OWN ROW IS MISSING!
    while read file
    do
        echo $file
        f2=${file##*/}
        CT=${f2%_*}
        sort +0 -1 $file | join $CIDIR/template_${mark} - -a 1 -e 0 | awk -v nam=$CT 'BEGIN{print nam}{out=($2>0?$2:0); print out}' | pr -mts $FULL - > $FULL.tmp

        # sort $file | awk -v nam=$CT 'BEGIN{print nam}{print $2}' | pr -mts $FULL - > $FULL.tmp
        mv $FULL.tmp $FULL
    done < <( ls ${DISTANCE_DIR}/BSS*_${mark}.txt | sort )
done < ${MARKS_LIST}

# Compress for sending and plotting:
cd ${DISTANCE_DIR}; tar -cvzf ${CIDIR}/distmat.tar.gz distance_all_*.tsv 

# Repeat for imputation directory:
cd ${IMPDIST_DIR}

rm ${IMPDIST_DIR}/imp_distance_all_*.tsv

# MERGE:
while read mark
do
    echo $mark
    FULL=${IMPDIST_DIR}/imp_distance_all_${mark}.tsv

    # List available:
    ls ${IMPDIST_DIR}/BSS*_${mark}.txt | awk -v mark=$mark '{sub(".*/","",$0); sub("_[A-Za-z0-9-]*.txt","",$0); print $0}' | sort +0 -1 > $CIDIR/template_${mark}

    # REMEMBER OWN ROW IS MISSING!
    while read file
    do
        echo $file
        f2=${file##*/}
        CT=${f2%_*}
        sort +0 -1 $file | join $CIDIR/template_${mark} - -a 1 -e 0 | awk -v nam=$CT 'BEGIN{print nam}{out=($2>0?$2:0); print out}' | pr -mts $FULL - > $FULL.tmp

        # sort $file | awk -v nam=$CT 'BEGIN{print nam}{print $2}' | pr -mts $FULL - > $FULL.tmp
        mv $FULL.tmp $FULL
    done < <( ls ${IMPDIST_DIR}/BSS*_${mark}.txt | sort )
done < ${MARKS_LIST}

# Compress for sending and plotting:
cd ${IMPDIST_DIR}; tar -cvzf ${CIDIR}/impdistmat.tar.gz imp_distance_all_*.tsv 


# Collate marks dist:
cd $MARKDIST_DIR

ALLMDIST=$CIDIR/mark_withinsamp_dist.tsv
rm $ALLMDIST
while read file; do
    fname=${file##*/}
    basename=${fname%%.txt}
    echo $basename
    mark=${basename%%_obs*}
    samp=${basename##*_}
    awk -v mark=$mark -v samp=$samp -vOFS="\t" '{print samp,mark, $1,$2}' $file >> $ALLMDIST
done < <(ls $MARKDIST_DIR/*.txt )

gzip -f $ALLMDIST
