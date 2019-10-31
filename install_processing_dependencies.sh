#!/bin/bash
# Install processing dependencies:
cd $SFTDIR

# SPP: 
if [[ "1" == "0" ]] 
then
    wget https://github.com/kundajelab/phantompeakqualtools/archive/master.zip -O phantompeakqualtools.zip 
    unzip phantompeakqualtools.zip

    Rscript -e '
    install.packages(c("caTools", "snowfall"), repos=:http://cloud.r-project.org/")
    source("http://bioconductor.org/biocLite.R")
    biocLite("Rsamtools")'

    # NOTE: make sure that my bzip2 is running (-fPIC compiled)
    # Drop .bzip2-1.0.6 from dotkits
    R CMD INSTALL phantompeakqualtools-master/spp_1.14.tar.gz
fi

# ==========
# Utilities:
# ==========
# BWtoBG:
if [[ ! -s $SFTDIR/bigWigToBedGraph ]]
then
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
    chmod +x bigWigToBedGraph
    mv bigWigToBedGraph $SFTDIR
fi

# BGtoBW:
if [[ ! -s $SFTDIR/bedGraphToBigWig ]]
then
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
    chmod +x bedGraphToBigWig
    mv bedGraphToBigWig $SFTDIR
fi

if [[ ! -s $SFTDIR/wigToBigWig ]]
then
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
    chmod +x wigToBigWig
    mv wigToBigWig $SFTDIR
fi


# bedtoBigBed:
if [[ ! -s $SFTDIR/bedToBigBed ]]
then
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
    chmod +x bedToBigBed
    mv bedToBigBed $SFTDIR
fi

# liftOver: 
if [[ ! -s $SFTDIR/liftOver ]]
then
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
    chmod +x liftOver
    mv liftOver $SFTDIR
fi

# bedClip: 
if [[ ! -s $SFTDIR/bedClip ]]
then
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedClip
    chmod +x bedClip
    mv bedClip $SFTDIR
fi


