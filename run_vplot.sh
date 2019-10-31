#!/bin/bash
# Grid Engine options
#$ -cwd
#$ -P compbio_lab
#$ -l h_vmem=8G 
#$ -l h_rt=1:00:00
#$ -tc 250 
#$ -j y
#$ -b y 
#$ -V 
#$ -r y 
#$ -o /broad/compbio/cboix/EPIMAP_ANALYSIS/db/out/ChIP
#$ -e /broad/compbio/cboix/EPIMAP_ANALYSIS/db/out/ChIP
# -------------------------------------
# Plot vplot of bam centered at bedfile
# -------------------------------------
# Load in variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

# NOTE: input ideally is PE; but we can try SE.
# Defaults:
EXTEND=1000
ALIGN="center"
NCORES=1
STRANDCOL=4
PRINTU=1
PRINTI=1
PRINTV=1
PRINTX=1
WINDOW=20
READLEN=-99

if [ $# -eq 0 ]; then
    echo >&2 "USAGE: qsub -t [RANGE] $0 [OPTIONS] 
    -i     Info table (tab-delimited) [required]
            - Input BAM file
            - Input BED file 
            - Output location
    -e     Number of bases to extend to each side [default: $EXTEND]
    -p     Alignment: center or ends (for paired end) [default: $ALIGN]
    -c     Number of threads to use [default: $NCORES]
    -s     Column # for BED file strand information [default: $STRANDCOL]
    -l     Set fixed template length (negative is ignored) [default: $READLEN]
    -u     Print uncompressed output file (0/1) [default: $PRINTU]
    -v     Print profile around bed file (0/1) [default: $PLOTU]
    -i     Print insert sizes across intervals (0/1) [default: $PRINTI]
    -x     Print sorted examples of top intervals (0/1) [default: $PRINTX]
    -w     Window size for plotting [default: $WINDOW]"
    exit 1
fi

while getopts i:e:p:c:s:u:v:i:w: o
do      case "$o" in
    i)		INFOTABLE="$OPTARG";;
    e)		EXTEND="$OPTARG";;
    p)		ALIGN="$OPTARG";;
    c)		NCORES="$OPTARG";;
    s)		STRANDCOL="$OPTARG";;
    l)      READLEN="$OPTARG";;
    u)		PRINTU="$OPTARG";;
    v)		PRINTV="$OPTARG";;
    i)		PRINTI="$OPTARG";;
    x)		PRINTX="$OPTARG";;
    w)		WINDOW="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

# Get files from info table:
BAMFILE=$( sed "${SGE_TASK_ID}q;d" ${INFOTABLE} | awk -v FS="\t" '{print $1}' )
BEDFILE=$( sed "${SGE_TASK_ID}q;d" ${INFOTABLE} | awk -v FS="\t" '{print $2}' )
OUTFILE=$( sed "${SGE_TASK_ID}q;d" ${INFOTABLE} | awk -v FS="\t" '{print $3}' )

# Add printing arguments:
PRINTOPTS=""
if [[ "$PRINTU" == "1" ]];then
    PRINTOPTS="${PRINTOPTS} -u"
fi

if [[ "$PLOTU" == "1" ]];then
    PRINTOPTS="${PRINTOPTS} -v"
fi

if [[ "$PRINTI" == "1" ]];then
    PRINTOPTS="${PRINTOPTS} -i"
fi

if [[ "$PRINTX" == "1" ]];then
    PRINTOPTS="${PRINTOPTS} -x"
fi

TMP_DIR=${TMP}/vplot_regions_${RANDOM}
mkdir -p ${TMP_DIR}

# Get basename/turn into bam file if a tagAlign file.
filebasename=${BAMFILE##*/}
basename=${filebasename%.gz}
suffix=${basename##*.}
if [[ "$suffix" == "tagAlign" ]]; then
    TMPFILE=${TMP_DIR}/${basename%.tagAlign}.tmp.bam
    INPUTFILE=${TMP_DIR}/${basename%.tagAlign}.bam
    # sorted tagAlign to bam and fix TLEN in column 9:
    bedtools bedtobam -i ${BAMFILE} -g ${CHROMSIZES} > ${TMPFILE} 
    samtools view -h $TMPFILE | awk -vOFS="\t" '{$9=$6; sub("M","", $9); print $0}' | samtools view -b -h | samtools sort -T ${TMP_DIR}/tmp_sort > ${INPUTFILE}
    samtools index ${INPUTFILE}
    repl="none"
else
    # Ensure sorted:
    INPUTFILE=${TMP_DIR}/${basename}
    samtools sort -T ${TMP_DIR}/tmp_sort -o ${INPUTFILE} ${BAMFILE}
    samtools index ${INPUTFILE}
    repl="chr"
fi

source activate mv_env
cmd="python ${BINDIR}/plot_vplot.py -a ${INPUTFILE} -b ${BEDFILE} -e $EXTEND -c $NCORES -o ${OUTFILE} -p ${ALIGN} ${PRINTOPTS} -s $STRANDCOL --window $WINDOW -r $repl -l $READLEN"
echo "[STATUS] Running command:"
echo "$cmd"
bash -c "$cmd"
source deactivate

# Clean up directories:
if [[ -d ${TMP_DIR} ]]; then 
    rm -rf ${TMP_DIR}
fi

end=`date +%s`
runtime=$((end-start))
echo "Finished vplot plotting and printing in $runtime seconds."

# === Test 1 (EZH2 in K562 to pks from K562)
# BAMFILE=/seq/epiprod/Data/Alignment_Post_Processing/007/73/Alignment_Post_Processing_773.bam
# BAMFILE=/broad/compbio/cboix/CRs/db/ChIP-seq/files/EZH2/tagAlign/FINAL_EZH2_BSS00264.sub.tagAlign.gz
# BEDFILE=/broad/compbio/cboix/CRs/db/ChIP-seq/files/EZH2/peaks/FINAL_EZH2_BSS00264.sub_VS_FINAL_WCE_BSS00264.narrowPeak.gz
# OUTFILE=/broad/compbio/cboix/CRs/db/vplot_EZH2_K562_vs_K562_pks
# STRANDCOL=6

# === Test 2 (EZH2 in K562 to pks from H1)
# BAMFILE=/seq/epiprod/Data/Alignment_Post_Processing/007/73/Alignment_Post_Processing_773.bam
# BAMFILE=/broad/compbio/cboix/CRs/db/ChIP-seq/files/EZH2/tagAlign/FINAL_EZH2_BSS00264.sub.tagAlign.gz
# BEDFILE=/broad/compbio/cboix/CRs/db/ChIP-seq/files/EZH2/peaks/FINAL_EZH2_BSS00177.sub_VS_FINAL_WCE_BSS00177.narrowPeak.gz
# OUTFILE=/broad/compbio/cboix/CRs/db/vplot_EZH2_K562_vs_H1_pks
# STRANDCOL=6
