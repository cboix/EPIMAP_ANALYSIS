#!/bin/bash
# Get and process mappability tracks (NOTE: requires processing deps)
cd $ANNDIR
if [[ ! -s $MAP36 ]] 
then
    wget https://www.pmgenomics.ca/hoffmanlab/proj/bismap/raw/hg19.umap.tar.gz
    tar -xvzf hg19

    gunzip hg19/k36.umap.wg.gz
    qsub -l h_rt=03:00:00 -l h_vmem=60G -N bwUMAPtrack -cwd -j y -b y -V -r y -o $DBDIR/out "wigToBigWig hg19/k36.umap.wg $CHROMSIZES k36.umap.bw; 
    bigWigToBedGraph k36.umap.bw k36.umap.bedGraph;
    awk 'BEGIN {FS=OFS="\t"} {if ($4 >= .75) print $1, $2, $3}' k36.umap.bedGraph | bedtools merge -i - > $MAP36;
    rm k36.umap.bw k36.umap.bedGraph"

    gunzip hg19/k50.umap.wg.gz
    qsub -l h_rt=03:00:00 -l h_vmem=60G -N bwUMAPtrack -cwd -j y -b y -V -r y -o $DBDIR/out "wigToBigWig hg19/k50.umap.wg $CHROMSIZES k50.umap.bw; 
    bigWigToBedGraph k50.umap.bw k50.umap.bedGraph;
    awk 'BEGIN {FS=OFS="\t"} {if ($4 >= .75) print $1, $2, $3}' k50.umap.bedGraph | bedtools merge -i - > $MAP50;
    rm k50.umap.bw k50.umap.bedGraph"

    gunzip hg19/k100.umap.wg.gz
    qsub -l h_rt=03:00:00 -l h_vmem=60G -N bwUMAPtrack -cwd -j y -b y -V -r y -o $DBDIR/out "wigToBigWig hg19/k100.umap.wg $CHROMSIZES k100.umap.bw; 
    bigWigToBedGraph k100.umap.bw k100.umap.bedGraph;
    awk 'BEGIN {FS=OFS="\t"} {if ($4 >= .75) print $1, $2, $3}' k100.umap.bedGraph | bedtools merge -i - > $MAP100;
    rm k100.umap.bw k100.umap.bedGraph"
fi

if [[ ! -s $BLACKLIST ]]; then
    wget https://www.encodeproject.org/files/ENCFF000KJP/@@download/ENCFF000KJP.bigBed -O $BLACKLIST
fi
