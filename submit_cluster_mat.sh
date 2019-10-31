#!/bin/bash
# ===================================
# Submit kmeans clustering for 
# binary matrices, in particular for
# prom/enh/dnase matrices
# ===================================
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh
export GWASJAR=${BINDIR}/gwas_enrichment/StateMWGwasPeakHyper.jar

# Enhancers:
K=300
MODELS=( "n15" "n18" )
TYPES=( "enh" "prom" )
model="n15"
fset="prom"

for model in ${MODELS[@]}; do
    for fset in ${TYPES[@]}; do
        DATADIR=${CHMMDIR}/calls/$model/$fset
        mkdir -p ${DATADIR}/clust/
        qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=5:00:00 -N cluster_${model}_${fset}_${K} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate pytorch_env > /dev/null; python ${BINDIR}/cluster_binary_mat.py -f ${DATADIR}/${model}_${fset}_bin_allchr_csr.cp -o ${DATADIR}/clust/autocls_${K} -K ${K}"
    done
done

# Main marks:
while read mark; do
    DATADIR=${CHMM_FMTDIR}/${mark}_binary
    mkdir -p ${DATADIR} ${DATADIR}/clust/
    qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=5:00:00 -N cluster_${mark}_${K} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate pytorch_env > /dev/null; python ${BINDIR}/cluster_binary_mat.py -f ${DATADIR}/${mark}_all_bin_allchr_csr -o ${DATADIR}/clust/autocls_${K} -K ${K}"
done < <( echo "DNase-seq\nH3K4me1\nH3K27ac")

# Marks:
while read mark; do
    DATADIR=${CHMM_FMTDIR}/${mark}_binary
    mkdir -p ${DATADIR} ${DATADIR}/clust/
    mv ${CHMM_FMTDIR}/${mark}_all_* ${DATADIR}
    # qsub -cwd -P compbio_lab -l h_vmem=20G -l h_rt=5:00:00 -N cluster_${mark}_${K} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate pytorch_env > /dev/null; python ${BINDIR}/cluster_binary_mat.py -f ${DATADIR}/${mark}_all_bin_allchr_csr -o ${DATADIR}/clust/autocls_${K} -K ${K}"
done < $MARKS_LIST


# ------------------
# MOTIF ENRICHMENT: 
# PROMOTER DNASE ENHANCER, 
# ALL AT K=300
# ------------------
# ENHANCERS
export MVDIR=${MAINDIR}/MOTIF_VALIDATION/
K=300
model="n15"
fset="enh"
prefix=cls_${K}
fpref=${model}_${fset}_k${K}_
DATADIR=${CHMMDIR}/calls/$model/$fset
CLSDIR=${DATADIR}/clust
DPREFIX=$CLSDIR/${fpref}${prefix}/${fpref}
OUTDIR=${CLSDIR}/gwas_enrichment
OUTPREFIX=${OUTDIR}/${fpref}
BEDFILE=${CLSDIR}/${prefix}_assignments.bed
mkdir -p ${CLSDIR}/${fpref}${prefix}/ $OUTDIR

# Fix bedfile:
# NOTE: Have to set to 200 because regions are of 200bp (cuts off end of genome).
sort -k1V ${CHROMSIZES_noY} | awk -vOFS="\t" 'BEGIN{a=0}; {print $1, a; n200 = 200 * int($2 / 200); a=a+n200}' | join $BEDFILE - | awk -vOFS="\t" '{print $1,$2-$5,$3-$5,$4}' > ${BEDFILE%bed}fixed.bed

# Make file directory:
echo "[STATUS] Separating files for prefix: ${fpref}"
awk '$1 == "chr17"' ${BEDFILE%bed}fixed.bed | head -n 2
awk -vOFS="\t" -v prefix=${CLSDIR}/${fpref}$prefix -v fpref=${fpref} '{print $1,$2,$3 > prefix"/"fpref$4".bed"}' ${BEDFILE%bed}fixed.bed

# Using Kmeans output, submit motif pipeline:
$MVDIR/submit_enrichment_motifs.sh -i ${CLSDIR}/${fpref}${prefix}/ -o "hg19"

# Run GWAS command:
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=3:00:00 -N gwas_${fpref} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${GWASJAR} $DPREFIX $OUTPREFIX $GWASCATALOG"





# PROMOTERS:
export MVDIR=${MAINDIR}/MOTIF_VALIDATION/
K=300
model="n15"
fset="prom"
prefix=autocls_${K}
fpref=${model}_${fset}_k${K}_
DATADIR=${CHMMDIR}/calls/$model/$fset
CLSDIR=${DATADIR}/clust
DPREFIX=$CLSDIR/${fpref}${prefix}/${fpref}
OUTDIR=${CLSDIR}/gwas_enrichment
OUTPREFIX=${OUTDIR}/${fpref}
BEDFILE=${CLSDIR}/${prefix}_assignments.bed
mkdir -p ${CLSDIR}/${fpref}${prefix}/ $OUTDIR

# Fix bedfile:
sort -k1V ${CHROMSIZES_noY} | awk -vOFS="\t" 'BEGIN{a=0}; {print $1, a; n200 = 200 * int($2 / 200); a=a+n200}' | join $BEDFILE - | awk -vOFS="\t" '{print $1,$2-$5,$3-$5,$4}' > ${BEDFILE%bed}fixed.bed

# Make file directory:
echo "[STATUS] Separating files for prefix: ${fpref}"
awk '$1 == "chr17"' ${BEDFILE%bed}fixed.bed | head -n 2
awk -vOFS="\t" -v prefix=${CLSDIR}/${fpref}$prefix -v fpref=${fpref} '{print $1,$2,$3 > prefix"/"fpref$4".bed"}' ${BEDFILE%bed}fixed.bed

# Using Kmeans output, submit motif pipeline:
$MVDIR/submit_enrichment_motifs.sh -i ${CLSDIR}/${fpref}${prefix}/ -o "hg19"

# Run command:
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=3:00:00 -N gwas_${fpref} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${GWASJAR} $DPREFIX $OUTPREFIX $GWASCATALOG"





# FOR DNASE:
# Separate out into chunks:
export MVDIR=${MAINDIR}/MOTIF_VALIDATION/
K=300
mark="DNase-seq"
prefix=autocls_${K}
fpref=${mark}_k${K}_
DATADIR=${CHMM_FMTDIR}/${mark}_binary
CLSDIR=${DATADIR}/clust
# Needs prefix and outprefix
DPREFIX=$CLSDIR/${fpref}${prefix}/${fpref}
OUTDIR=${CLSDIR}/gwas_enrichment
OUTPREFIX=${OUTDIR}/${fpref}
BEDFILE=${CLSDIR}/${prefix}_assignments.bed
mkdir -p ${CLSDIR}/${fpref}${prefix}/ ${OUTDIR}/ 

# Fix bedfile:
sort -k1V ${CHROMSIZES_noY} | awk -vOFS="\t" 'BEGIN{a=0}; {print $1, a; n200 = 200 * int($2 / 200); a=a+n200}' | join $BEDFILE - | awk -vOFS="\t" '{print $1,$2-$5,$3-$5,$4}' > ${BEDFILE%bed}fixed.bed

# Make file directory:
echo "[STATUS] Separating files for prefix: ${fpref}"
awk '$1 == "chr17"' ${BEDFILE%bed}fixed.bed | head -n 2
awk -vOFS="\t" -v prefix=${CLSDIR}/${fpref}$prefix -v fpref=${fpref} '{print $1,$2,$3 > prefix"/"fpref$4".bed"}' ${BEDFILE%bed}fixed.bed

# Using Kmeans output, submit motif pipeline:
$MVDIR/submit_enrichment_motifs.sh -i ${CLSDIR}/${fpref}${prefix}/ -o "hg19"

# Run command:
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=3:00:00 -N gwas_${fpref} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${GWASJAR} $DPREFIX $OUTPREFIX $GWASCATALOG"

