#!/bin/bash
# Get ROADMAP Data
export CHMM=$DBDIR/Roadmap/ChromHMM15state
mkdir -p $CHMM 

if [[ ! -s $CHMM/E001_15_coreMarks_dense.bed.gz ]] 
then
    cd $CHMM
    wget http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.dense.browserFiles.tgz
    tar -xvzf all.dense.browserFiles.tgz
    rm all.dense.browserFiles.tgz
fi

# Get enhancer parts of the track: 
while read denseBED
do
    echo "Processing file $denseBED"
    enhBED=$(echo ${denseBED} | sed -r -e 's/\dense\.bed/enhancers\.bed/g')
    if [[ ! -s $enhBED ]] 
    then
        zcat $denseBED | awk '$4 ~ /Enh/{print $0}' | gzip -c > $enhBED
    fi
done < <( ls ${CHMM}/*15_coreMarks_dense.bed.gz )

