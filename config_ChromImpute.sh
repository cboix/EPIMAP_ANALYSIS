#!/bin/bash
# Configuration file for running ChromImpute and downstream analyses
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]; then
    export MAINDIR=/broad/compbio/cboix
    export SFTDIR=${MAINDIR}/software
else
    export MAINDIR=$HOME
    export SFTDIR=${HOME}/local
fi
export DBDIR=${MAINDIR}/EPIMAP_ANALYSIS/db
export BINDIR=${MAINDIR}/EPIMAP_ANALYSIS/bin
export IMGDIR=${MAINDIR}/EPIMAP_ANALYSIS/img
export ANNDIR=${DBDIR}/Annotation
export SNDIR=${BINDIR}/snake_pipe
export OUTDIR=$DBDIR/out/ChromImpute
mkdir -p ${DBDIR}/out $OUTDIR $ANNDIR ${DBDIR}/out/annotate
start=`date +%s`
hostname -f

if [[ $# -lt 1 ]];then
    NUMSTATES=18
    echo "NUMSTATES is ${NUMSTATES} (default)." 
else 
    NUMSTATES=$1
    echo "Set NUMSTATES to ${NUMSTATES}." 
fi

# config Variables:
export MEM=8  # Memory for ChromHMM runs
export LD_LIBRARY_PATH=${SFTDIR}/bzip2-1.0.6:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${SFTDIR}/miniconda2/envs/mv_env/lib/:${LD_LIBRARY_PATH}

# ChromImpute and ChromHMM:
# NOTE: If malfunctioning, remake using: make clean; make
export CHROMIMPUTE=${BINDIR}/ChromImpute_SRC/ChromImpute.jar
export CHROMHMM=${SFTDIR}/ChromHMM/ChromHMM.jar
export GWASJAR=${BINDIR}/gwas_enrichment/StateMWGwasPeakHyper.jar
export MVDIR=${MAINDIR}/MOTIF_VALIDATION/

# Genome Annotation Files (all for hg19)
export MAPDIR=$ANNDIR/umap
export GTFDIR=$ANNDIR/GENCODE
export CHROMSIZES_noY=$ANNDIR/hg19.chrom.sizes_noY
export CHROMSIZES=$ANNDIR/hg19.chrom.sizes_main
export CHROMARM=$ANNDIR/chromArm.bed
export LOCHAIN=$ANNDIR/hg38ToHg19.over.chain.gz
export GENCODE=$GTFDIR/gencode.v27lift37.primary_assembly.annotation.gtf.gz
export GTF_PREF=${GENCODE#*.}
export GTF_SUFFIX=${GTF_PREF%%.annotation.gtf.gz}
export GENE_GTF=$GTFDIR/Gene.${GTF_SUFFIX}.bed
export GENE_COLS=$GTFDIR/Gene.${GTF_SUFFIX}.cols.tsv
mkdir -p $GTFDIR $MAPDIR

# For GWAS/LDAK/LDSC:
export GENOMEDIR=${MAINDIR}/genomes
export LDSCBASE=${GENOMEDIR}/ldsc_baseline
# export GWRECENTCAT=$DBDIR/gwas_catalog_v1.0-associations_e93_r2018-12-21.tsv
# export GWASCATALOG=${DBDIR}/gwascatalog_dec21_2018.txt
export GWRECENTCAT=$DBDIR/gwas_catalog_v1.0-associations_e96_r2019-05-03.tsv 
export GWASCATALOG=${DBDIR}/gwascatalog_may03_2019.txt
export H3SNPS=${GENOMEDIR}/hapmap3/w_hm3.snplist
export LDSCDIR=$SFTDIR/ldsc
if [[ ! -s $GWASCATALOG ]];then
    if [[ ! -s $GWRECENTCAT ]]; then
        # Get 20190503 annotation
        wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/2019/05/03/gwas-catalog-associations.tsv -O $DBDIR/gwas_catalog_v1.0-associations_e96_r2019-05-03.tsv
    fi
    zcat ${DBDIR}/gwascatalog_sep12_2014.txt.gz | head -n 1 > $GWASCATALOG
    awk -F"\t" -vOFS="\t" 'NR>1 && $12 !~/ x /{split($12,a,";"); print "1", a[1], $13, $13 + 1, $22, $2, $3, $4, $5, $7, $8, $9, $10, $11, $15, $21, $27, $28, $30, $31, $32, $33, $34}' $GWRECENTCAT >> $GWASCATALOG
    # TAKING FIRST LOCATION IN STRINGS OF MANY
    awk -F"\t" -vOFS="\t" '{gsub("\t\t", "\tNR\t",$0); print $0}' $GWASCATALOG | awk -F"\t" -vOFS="\t" '{gsub("[-_;].*","",$3); gsub("[-_;].*","",$4); gsub("[-_;].*","",$5); gsub("[-_;].*","",$15); gsub("[-_;].*","",$16); print $0}' > $GWASCATALOG.tmp
    # Fix about ~3k that have NR (not reported - in $5) - fixed about 2k.
    awk -F"\t" -vOFS="\t" '{if($2 == "NR"){chr=split($5,a,":"); sub("chr","",a[1]); sub("Chr","",a[1]); if(a[1] ~ /^[0-9XY][0-9XY]?$/){$2=a[1]; $3=a[2];}}; if($2 == "23"){$2="X"}; $4=$3 + 1; print $0}' $GWASCATALOG.tmp | awk -F"\t" '$2 != "NR"{print $0}' > $GWASCATALOG
    # Sanity checks:
    awk -F"\t" -vOFS="\t" '$2 == "NR"{print $5}' $GWASCATALOG | sort
    awk -F"\t" -vOFS="\t" '{print $2}' $GWASCATALOG | sort | uniq -c | sort -nr
    # LiftOver coordinates to hg19 - lose about 0.1% that had errors + other issues.
    awk -F"\t" -vOFS="\t" 'NR > 1{print "chr"$2,$3,$4, NR - 1}' $GWASCATALOG > ${GWASCATALOG}_coords.bed
    liftOver ${GWASCATALOG}_coords.bed ${LOCHAIN} ${GWASCATALOG}_hg19_coords.bed ${GWASCATALOG}_coords.unmapped
    # Join back with original (all that can be matched):
    awk -vOFS="\t" '{print NR - 1, $0}' $GWASCATALOG | sort +0 -1 > $GWASCATALOG.tmp
    awk -vOFS="\t" 'BEGIN{print 0,"chrom", "chromStart"}{sub("chr","",$1); print $4,$1,$2}' ${GWASCATALOG}_hg19_coords.bed | sort +0 -1 | join - $GWASCATALOG.tmp -t $'\t' | awk -F"\t" -vOFS="\t" '{gsub("\%","pct",$0); $5=$2; $6=$3; $7=$6+1; if(NR == 1){$7="chromEnd"}; printf $4; for(i=5; i<=NF; i++){a=($i==""?"NR":$i); printf "\t"a}; printf "\n"}' > $GWASCATALOG
    rm $GWASCATALOG.tmp
fi

# ChromImpute directories:
export CIDIR=$DBDIR/ChromImpute
export CONVERTED_DIR=$CIDIR/converted
# Distance directories:
export DISTANCE_DIR=$CIDIR/distance
export REGDIST_DIR=$CIDIR/region_distance
export IMPDIST_DIR=$CIDIR/imputed_distance
export MIXDIST_DIR=$CIDIR/miximpobs_distance
export DIFFDIST_DIR=$CIDIR/difference_distance
export MARKDIST_DIR=$CIDIR/mark_distance
export DISTR_DIR=$CIDIR/signal_distribution
# Other directories
export TRAIN_DIR=$CIDIR/train_data
export PREDICTOR_DIR=$CIDIR/predictor
export IMPUTED_DIR=$CIDIR/imputed
export IMPSAMPLE_DIR=$CIDIR/imp_xsample
export PKIMP_DIR=$CIDIR/imputed_peaks
export PKOBS_DIR=$CIDIR/observed_peaks
export MPD_DIR=$CIDIR/mergedpeak_data
export EVAL_DIR=$CIDIR/stats
export CHECK_DIR=$CIDIR/checks
export IMPFMT_DIR=$CIDIR/binarized_imputed
export IMPFMTQCUT_DIR=$CIDIR/binarized_imputed_qcutoff
export SUBSET_DIR=$CIDIR/subsetted_tracks
mkdir -p ${CIDIR} ${CONVERTED_DIR} ${DISTANCE_DIR} ${TRAIN_DIR} ${CHECK_DIR}
mkdir -p ${PREDICTOR_DIR} ${IMPUTED_DIR} ${EVAL_DIR} ${IMPFMT_DIR} ${IMPFMTQCUT_DIR} ${IMPSAMPLE_DIR}
mkdir -p $MIXDIST_DIR $REGDIST_DIR ${IMPDIST_DIR} ${DIFFDIST_DIR} ${PKIMP_DIR} ${PKOBS_DIR} ${MPD_DIR} ${MARKDIST_DIR} ${DISTR_DIR} ${SUBSET_DIR}

# Files:
export MAPFILE=${ANNDIR}/all_submitted_released_biosample_mapping.tsv
export KEPTCELLS=${ANNDIR}/kept_bssid_20190322.txt

# ChromHMM directories
export AMODELDIR="/broad/compbio/anshul/projects/roadmap/chromhmmSegmentations/ChmmModels"
export CHMMDIR=$DBDIR/ChromHMM
export RDCHMMDIR=${DBDIR}/roadmap_chromHMM
export CHMM_FMTDIR=${CHMMDIR}/binarized
export BYCELLDIR=${CHMMDIR}/files_by_cell
export CALLDIR=${CHMMDIR}/calls
export FREQDIR=${CHMMDIR}/freqs
export EXPTINFO=${CHMMDIR}/expt_list.tsv
export CELLINFO=${CHMMDIR}/cell_list.tsv
mkdir -p ${CHMMDIR} ${CHMM_FMTDIR} ${FILEDIR} ${CALLDIR} $RDCHMMDIR

# List all CHMM experimental data:
export DATATAG="c_t.sub"
if [[ ! -s $EXPTINFO ]] || [[ ! -s $CELLINFO ]]; then
    cd ${CHMM_FMTDIR}
    ls */BSS[0-9]*_${DATATAG}_chr1_binary.txt.gz | awk -vOFS="\t" '{a=$0; sub("/","\t",a); sub("_c.*gz","",a); sub("_chr1.*gz","",$0); print a,$0}' | awk -vOFS="\t" '{print $2,$1,$3}' | sort -u > $EXPTINFO
    # NOTE: There are about 40 cells with CTCF only. They will be fully imputed.
    awk '$2 ~ /^H[1234]/||/CTCF/||/DNase-seq/||/ATAC-seq/{print $1}' ${EXPTINFO} | sort -u > ${CELLINFO}
    cd $CHMMDIR
fi

# Public directories:
export PUBLIC_BWDIR=$DBDIR/public_bigwigs_released
export PUBLIC_BDGDIR=$DBDIR/public_bedgraphs_released
export PUBLIC_CMDIR=$DBDIR/public_ChromHMM_released
export PUBLIC_METADIR=$DBDIR/public_metadata_released
mkdir -p $PUBLIC_BWDIR $PUBLIC_CMDIR $PUBLIC_BWDIR/observed/ $PUBLIC_BWDIR/imputed/ $PUBLIC_METADIR $PUBLIC_BDGDIR $PUBLIC_BDGDIR/observed/ $PUBLIC_BDGDIR/imputed/

# Accessory epitopes: 
accs_epitopes=(CTCF EP300 POLR2A SMC3 RAD21)
main_epitopes=(DNase-seq H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
imp_epitopes=(DNase-seq H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 H2AFZ H3K4me2 H3K79me2 H3K9ac H4K20me1)

# Get ChromHMM model:
if [[ "$NUMSTATES" == "18" ]]; then 
    MODEL_LOC="${AMODELDIR}/core_K27ac/parallel/set1/final/model_18_core_K27ac.txt"
    MODELNAME="observed_aux_18";
    #epitopes=(H3K4me1 H3K4me3 H3K27me3 H3K36me3 H3K9me3 H3K27ac)
    # Order in the model file:
    epitopes=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
elif [[ "$NUMSTATES" == "15" ]]; then
    MODEL_LOC="${AMODELDIR}/coreMarks/parallel/set2/final/model_15_coreMarks.txt"
    MODELNAME="observed_15";
    # epitopes=(H3K4me3 H3K4me1 H3K27me3 H3K9me3 H3K36me3)
    # Order in the model file:
    epitopes=(H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
elif [[ "$NUMSTATES" == "25" ]]; then
    # NOTE: Real (non-symlinked) 25 state model is in
    # /broad/compbio/jernst/compbio-hp/SIGNALPREDICT6/TIER2_REORDER_MODEL/states25
    MODEL_LOC="$AMODELDIR/imputed12marks/n25/model_25_imputed12marks.txt"
    MODELNAME="imputed12_25";
    epitopes=(DNase-seq H3K27ac H3K9ac H3K4me3 H3K4me1 H3K36me3 H3K9me3 H3K27me3 H3K79me2 H4K20me1 H3K4me2 H2AFZ)
else
    echo "Not recognized setting for $NUMSTATES"
    echo "Valid values are 15, 18, and 25"
fi

# Model-specific directories:
export NCALLDIR=${CALLDIR}/${MODELNAME}
export STATEDIR=${NCALLDIR}/STATEBYLINE
export MODEL_FILE=${ANNDIR}/ChromHMM_model_${MODELNAME}_states.txt
cat $MODEL_LOC > $MODEL_FILE
mkdir -p ${NCALLDIR} ${STATEDIR} ${FREQDIR}

# Information tables:
export SAMPLEMARK_TAB=${CIDIR}/sample_mark_table.tsv
export EXT_SAMPMARK_TAB=${SAMPLEMARK_TAB%.tsv}_extended.tsv
export TOIMPUTE_TAB=${CIDIR}/imputation_table.tsv
export ACCSIMPUTATION_TAB=${CIDIR}/accs_imputation_table.tsv
export EXTRAIMPUTATION_TAB=${CIDIR}/extra_imputation_table.tsv
export FULLIMPUTATION_TAB=${CIDIR}/full_imputation_table.tsv
export MARKS_LIST=${CIDIR}/marks_available.tsv
export IMPOBS_TAB=${CIDIR}/impobs_table_fordist.tsv
export ALL_TRACKS_TAB=${CIDIR}/all_impobs_tracks_table.tsv
export ALL_UQ_TAB=${CIDIR}/all_impobs_tracks_uniq_table.tsv
export MIXOBS_TAB=${CIDIR}/mixobs_table_fordist.tsv
export DIFF_TAB=${CIDIR}/obsimp_diff_table_fordist.tsv

# DNase masterlist:
export DML_DIR=${DBDIR}/DHS_Index_WM201902
export DMLPREF=${DML_DIR}/masterlist_DHSs_733samples_WM20180608_all_coords
export MLCOORD=${DMLPREF}.txt
export LOCOORD=${DMLPREF}_hg19.txt
export CORECOORD=${DMLPREF}_hg19.core.srt.txt
export NUMCOORD=${CORECOORD%txt}numbered.txt
mkdir -p ${DML_DIR}
# Ensure masterlist coords are in hg19:
if [[ ! -s $LOCOORD ]] || [[ ! -s $CORECOORD ]]; then
    liftOver ${MLCOORD} ${LOCHAIN} ${LOCOORD} $DML_DIR/ml.over.unmapped
    # Take core only + sort (separately)
    awk '$1 ~ /^chr[0-9X][0-9]?$/' $LOCOORD | sort -k1,1V > ${DMLPREF}_hg19.core.txt
    sort -k1,1V -k2n ${DMLPREF}_hg19.core.txt > ${CORECOORD}
    # Number 
    awk '{print NR, $1"-"$2"-"$3}' $CORECOORD | sort +0 -1 > ${NUMCOORD}
    # Print number of regions in each:
    wc -l $MLCOORD
    wc -l $LOCOORD
    wc -l $CORECOORD
fi

# Enrichment directories:
export GO_DIR=${DBDIR}/go_enrichments
export GWAS_DIR=${DBDIR}/gwas_enrichments
export MOTIF_DIR=${DBDIR}/Enrichments
mkdir -p ${GO_DIR} ${GWAS_DIR} ${MOTIF_DIR}

# Genomes for motifs:
export GENOME="b37"
export GENOME_PREF="human_g1k_v37"
export REF_GATK=${GENOMEDIR}/gatk_resource_bundle_${GENOME}
export REF_GENOME_FULL=${REF_GATK}/${GENOME_PREF}.fasta  # Primary + decoys + alt

if [[ ! -s $CHROMSIZES ]]; then
    wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/chromInfo.txt.gz -O $CHROMSIZES\.gz
    gunzip -c $CHROMSIZES.gz | awk -F"\t" -vOFS="\t" '$1 ~ /chr[0-9XY]+$/{print $1,$2}' > $CHROMSIZES
    rm $CHROMSIZES.gz
fi

# Run ID (if none):
if [[ "$UQID" == "" ]]; then
    export UQID=ChromImpute_$RANDOM
    export TODAY=$(date +%Y%m%d-%H%M)
    echo ${UQID}
fi
