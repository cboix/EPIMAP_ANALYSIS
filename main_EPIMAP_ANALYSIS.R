#!/bin/bash
# Get ROADMAP Data
# Directories:
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]
then
    export DBDIR=$HOME/data/EPIMAP_ANALYSIS/db
    export BINDIR=$HOME/data/EPIMAP_ANALYSIS/bin
    export SFTDIR=$HOME/data/software
else
    export DBDIR=$HOME/EPIMAP_ANALYSIS/db
    export BINDIR=$HOME/EPIMAP_ANALYSIS/bin
    export SFTDIR=$HOME/bin
fi
export ANNDIR=$DBDIR/Annotation
mkdir -p $DBDIR/out $ANNDIR

# Variables:
export RUNALL=false # Recreate everything?
export sleep_time=30 # Use if we need bash traps.

# For filtering reads (not necessary now)
# export FILTDIR=$HOME/cHMM/bin
# export SEQDIR="/broad/compbio/anshul/projects/encode/rawdata/sequence"
# export UMAPDIR="/broad/compbio/anshul/projects/umap" 
# export SDIR="${SEQDIR}/hg19"
# export UDIR="${UMAPDIR}/hg19/globalmap_k20tok54"

# Software + Directories:
export SAMTL="/broad/software/free/Linux/redhat_6_x86_64/pkgs/samtools/samtools_1.3.1/bin/samtools"
export SPPDIR="/broad/compbio/cboix/software/phantompeakqualtools-master"
export CHROMIMPUTE="/broad/compbio/cboix/software/ChromImpute/ChromImpute.jar"
export MAPDIR=$ANNDIR/umap
export GTFDIR=$ANNDIR/GENCODE
mkdir -p $GTFDIR $MAPDIR

# Data directories:
export VLDIR=$DBDIR/Validation
export SAMPDIR=$DBDIR/Validation/Sampling
export IHECDIR=$DBDIR/IHEC
export CHPDIR=${DBDIR}/ChIP-seq
export CHPLNK=${CHPDIR}/file_links/all_submitted_released
export DHSDIR=${DBDIR}/DNase-seq
export DHSINFO=${DHSDIR}/file_links/all_submitted_released/DNase-seq_bam.tsv
export FILE_DIR=$DHSDIR/files/tagAlign
export WINDIR=$DHSDIR/files/window_coverage
export QCDIR=$DHSDIR/files/qc
mkdir -p $FILE_DIR $QCDIR $WINDIR $IHECDIR $CHPDIR $VLDIR $SAMPDIR

# For ChromImpute:
export CIDIR=$DBDIR/ChromImpute
export CONVERTED_DATADIR=$CIDIR/converted
mkdir -p $CIDIR ${CONVERTED_DATADIR}

# Genome Annotation Files: 
export CHROMSIZES=$ANNDIR/hg19.chrom.sizes
export CHROMARM=$ANNDIR/chromArm.bed
export LOCHAIN=$ANNDIR/hg38ToHg19.over.chain.gz
export GENCODE=$GTFDIR/gencode.v27lift37.primary_assembly.annotation.gtf.gz
export GTF_SUFFIX=${${GENCODE#*.}%%.annotation.gtf.gz} # othw serial - fails on bash. 
export GENE_GTF=$GTFDIR/Gene.${GTF_SUFFIX}.bed
export GENE_COLS=$GTFDIR/Gene.${GTF_SUFFIX}.cols.tsv

# Mappability regions
export MAP36=$MAPDIR/k36.Umap.MultiTrackMappability.filtered.bed
export MAP50=$MAPDIR/k50.Umap.MultiTrackMappability.filtered.bed
export MAP100=$MAPDIR/k100.Umap.MultiTrackMappability.filtered.bed
export BLACKLIST=$ANNDIR/hg19_blacklist_ENCFF000KJP.bigBed

# Get necessary files: 
if [[ ! -s $CHROMSIZES ]] 
then
    # wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -O $CHROMSIZES
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -O $CHROMSIZES
fi

if [[ ! -s $LOCHAIN ]]
then
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O $LOCHAIN
fi

if [[ ! -s $GENCODE ]] 
then
    wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.annotation.gtf.gz -O $GENCODE
    zcat $GENCODE | awk '$3 ~/gene/' > ${GENE_GTF}

    # Make gene columns:
    awk -vFS="\t" -vOFS="\t" '{sub("gene_id \"","",$0); sub("\";","\t",$0); print $1,$9 }' ${GENE_GTF} > ${GENE_COLS}
fi

if [[ ! -s $CHROMARM ]]
then
    cd $ANNDIR
    wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz -O $ANNDIR/cytoBand.txt.gz
    zcat $ANNDIR/cytoBand.txt.gz | grep acen | awk '$4 ~ /p/{print $1,$3}' | sort +0 -1 > $ANNDIR/centromeres_hg19.txt

    sort +0 -1 $CHROMSIZES | join - $ANNDIR/centromeres_hg19.txt | awk -vOFS="\t" '{print $1,0,$3,$1"A"; print $1,$3,$2,$1"B"}' >  $CHROMARM
fi


# Install spp + utilities if not installed
source $BINDIR/install_processing_dependencies.sh
source $BINDIR/make_TSS_windows.sh # Windows around genes (for genome cov etc.)
source $BINDIR/get_mappability_tracks.sh 

cd $DBDIR
# TODO qsub
# Get links from ENCODE directory + make matrices of availability
if [[ ! -s ${ANNDIR}/"experiments_all_submitted_released.RData" ]] 
then
    R --slave -f $BINDIR/get_experiments_ENCODE.R

    # TODO logic here
    R --slave -f $BINDIR/get_available_ENCODE_links.R ChIP-seq bam
    R --slave -f $BINDIR/get_available_ENCODE_links.R DNase-seq bam
    R --slave -f $BINDIR/get_available_ENCODE_links.R WGBS bam
    R --slave -f $BINDIR/get_available_ENCODE_links.R RRBS bam # TODO fix for RRBS -- data not available?
fi

# Read info table of all cell types:
awk '{ }' $ANNDIR/all_submitted_released_cell_types.tsv

# prioritization for mark processing (first are histone mods, followed by other ChIP-seq data)
awk '{ }' $ANNDIR/all_submitted_released_ChIP_abundance.tsv


# Get data + process to low-level ENH:
source $BINDIR/get_roadmap_enh.sh
source $BINDIR/get_DNase_data.sh

# =======================
# DHS pipe (file by file)
# =======================
DEXP=$(wc -l $DHSINFO | awk '{print $1}')
echo $DEXP 

# Run DHS pipeline: -t 2-$DEXP
if [[ "0" = "1" ]]
then
    # ONLY run if REALLY NECESSARY to do all (skip 1: gives header)
    # TODO NOT RIGHT PARAMS:
    qsub -cwd -t 2-$DEXP -l h_rt=24:00:00 -l h_vmem=40G -N process_DHS -j y -b y -V -r y -o $DBDIR/out $BINDIR/run_DNase_pipeline.sh 

    # For reruns on long queue:
    RERUN_INFO=${FILE_DIR}/RERUN_DNase-seq_bam.tsv
    head -n 1 $DHSINFO > ${RERUN_INFO}
    for SID in `seq 2 $DEXP`
    do
        run=$( sed "${SID}q;d" $DHSINFO )
        eval $( echo ${run} | awk '{sub("tered align","tered_align",$0); a = $1; sub(".*download/","",a); sub("\\..*$","",a); printf("id=\"%s\"; link=\"%s\"; assembly=\"%s\"; cell=\"%s\"; processing=\"%s\";",a,$1,$3,$5,$4)}' )
        FTA=${FILE_DIR}/FINAL_${id}_${cell}.tagAlign
        if LC_ALL=C gzip -l ${FTA} | awk 'NR==2 {exit($2!=0)}'
        then
            echo $SID
            echo $run >> ${RERUN_INFO}
        fi
    done

    # NOTE will need to run last steps for 2-64 again! (bad umap)
    # Run on long queue - needs more time to download/process.
    REXP=$(wc -l ${RERUN_INFO} | awk '{print $1}')
    qsub -cwd -t 2-$REXP -q long -l h_vmem=50G -N reprocess_DHS -j y -b y -V -r y -o $DBDIR/out $BINDIR/run_DNase_pipeline.sh ${RERUN_INFO}

    # TMP: Removing bad files and re-running:
    cd $WINDIR
    while read file
    do
        if [[ "9" != "$( awk 'NR==2{print NF}' $file )" ]] 
        then
            rm $file
            id=${file%%_*}
            cell=${${file%.DHS_windows}#*_}
            echo "Re-running for ${id}: $cell"
            source $BINDIR/DNase_compute_window_occupancy.sh $id $cell
        fi
    done < <( ls *.DHS_windows )

    # MERGING DHS FILES: 
    # Only files that are fully complete: 
    MINF=${FILE_DIR}/MERGING_DNase-seq.tsv
    rm -f $MINF
    cd $WINDIR
    while read file
    do
        if [[ "9" == "$( awk 'NR==2{print NF}' $file )" ]] 
        then
            id=${file%%_*}
            cell=${${file%.DHS_windows}#*_}
            echo "$file\t$id\t$cell" >> $MINF
        fi
    done < <( ls *.DHS_windows )

    for ((i = 3; i <= 9; i++))
    do
        echo $i
        # SET UP FILE: 
        WIN=$(awk -v ind=$i 'NR==1{print $ind}' ${WINDIR}/ENCFF000SGF_B_cell.DHS_windows)
        OUTFILE=${WINDIR}/all_${WIN}_raw.tsv
        if [[ ! -s $OUTFILE ]] 
        then
            echo "STARTING to make $OUTFILE"
            IFS=$'\t'
            while read file id cell
            do
                echo "Add $id, $cell"
                awk -v ind=$i -v id=$id -v cell=$cell 'NR==1{print id"_"cell}NR>1{print $ind}' $file | pr -mts $OUTFILE - > ${OUTFILE}.tmp
                mv ${OUTFILE}.tmp ${OUTFILE}
            done < $MINF
            echo "FINISHED making $OUTFILE"
        fi
    done

    # TODO Update the arsaux package on the broad server so that it handles arbitrary covariates!!

fi

# ====================================
# ChIP-seq pipe (an epitope at a time)
# ====================================
# =================================================
# TODO make master file to decides how to run each:
# long array for the small epitopes, 
# submit each of the histone marks + high amnts in a parallel manner.
# =================================================
# Should be number of epitopes:
# CEXP=$(wc -l $CHPINFO | awk '{print $1}')
epitope=H3K27me3
# TODO add SCCA, MACS, fc, and log10pval to pipeline
# NOTE WCE and CONTROL DO THEY DIFFER?
# eplist=(BACH1 BRF2 SIN3A HDAC3 CTBP2)
# for epitope in ${eplist[@]}
while read epitope
do
    echo $epitope
    FNUM=$(wc -l ${CHPLNK}/${epitope}_bam.csv | awk '{print $1 - 1}')
    if [[ "$FNUM" -lt "15" ]]
    then 
        qsub -cwd -q long -l h_vmem=30G -N process_${epitope} -j y -b y -V -r y -o $DBDIR/out/ChIP "$BINDIR/run_ChIP_epitope.sh $epitope"
    else 
        # Alternatively: 
        source $BINDIR/run_ChIP_parallel.sh $epitope
    fi
done < $DBDIR/TF_names2


# ===============
# IHEC PIPELINE: 
# ===============
# Get IHEC DATA: 
awk -vFS="\t" '{sub("_[IH][n3].*","",$6); print $6}' $ANNDIR/IHEC_bw_histonemods.tsv | sort -u > $ANNDIR/IHEC_donors.tsv
FNUM=$(wc -l ${CHPLNK}/${epitope}_bam.csv | awk '{print $1 - 1}')

IEXP=$(wc -l ${ANNDIR}/IHEC_donors.tsv | awk '{print $1}')
qsub -cwd -t 1-$IEXP -l h_vmem=30G -N proc_IHECdata -j y -b y -V -r y -o $DBDIR/out $BINDIR/get_IHEC_histone_bw_tracks.sh





