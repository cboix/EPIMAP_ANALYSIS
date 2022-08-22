#!/bin/bash
# ===================================
# Submit aggregate marks and chromhmm
# ===================================
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

# Enhancers:
# Aggregating jobs:
fset="enh"
model="15"
states="6,7,12"
DATADIR="${CHMMDIR}/calls/n${model}/$fset"
STATEDIR="${CHMMDIR}/calls/n${model}/STATEBYLINE/"
MAINJOB="source activate pytorch_env > /dev/null; python ${BINDIR}/aggregate_binary_txt_to_npy.py --dir ${STATEDIR} --out ${DATADIR}/n${model}_${fset}_bin --suffix _statebyline.txt.gz --mid _${model}_CALLS_PER_LINE_ --state $states --chromhmm"
while read chr size; do
    qsub -cwd -P compbio_lab -l h_vmem=25G -l h_rt=5:00:00 -N binarize_${model} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$MAINJOB --chr $chr"
done
qsub -cwd -P compbio_lab -l h_vmem=25G -l h_rt=2:00:00 -N binarize_${model}_all -hold_jid binarize_${model} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$MAINJOB"

fset="enh"
model="18"
states="7,8,9,10,11,15"
DATADIR="${CHMMDIR}/calls/n${model}/$fset"
STATEDIR="${CHMMDIR}/calls/n${model}/STATEBYLINE/"
MAINJOB="source activate pytorch_env > /dev/null; python ${BINDIR}/aggregate_binary_txt_to_npy.py --dir ${STATEDIR} --out ${DATADIR}/n${model}_${fset}_bin --suffix _statebyline.txt.gz --mid _${model}_CALLS_PER_LINE_ --state $states --chromhmm"
while read chr size; do
    qsub -cwd -P compbio_lab -l h_vmem=25G -l h_rt=5:00:00 -N binarize_${model} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$MAINJOB --chr $chr"
done
qsub -cwd -P compbio_lab -l h_vmem=25G -l h_rt=2:00:00 -N binarize_${model}_all -hold_jid binarize_${model} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$MAINJOB"


# Promoters:
fset="prom"
model="15"
states="1,2,10"
DATADIR="${CHMMDIR}/calls/n${model}/$fset"
STATEDIR="${CHMMDIR}/calls/n${model}/STATEBYLINE/"
MAINJOB="source activate pytorch_env > /dev/null; python ${BINDIR}/aggregate_binary_txt_to_npy.py --dir ${STATEDIR} --out ${DATADIR}/n${model}_${fset}_bin --suffix _statebyline.txt.gz --mid _${model}_CALLS_PER_LINE_ --state $states --chromhmm"
while read chr size; do
    qsub -cwd -P compbio_lab -l h_vmem=25G -l h_rt=5:00:00 -N binarize_p${model} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$MAINJOB --chr $chr"
done
qsub -cwd -P compbio_lab -l h_vmem=25G -l h_rt=2:00:00 -N binarize_p${model}_all -hold_jid binarize_p${model} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$MAINJOB"


fset="prom"
model="18"
states="1,2,3,4,14"
DATADIR="${CHMMDIR}/calls/n${model}/$fset"
STATEDIR="${CHMMDIR}/calls/n${model}/STATEBYLINE/"
MAINJOB="source activate pytorch_env > /dev/null; python ${BINDIR}/aggregate_binary_txt_to_npy.py --dir ${STATEDIR} --out ${DATADIR}/n${model}_${fset}_bin --suffix _statebyline.txt.gz --mid _${model}_CALLS_PER_LINE_ --state $states --chromhmm"
while read chr size; do
    qsub -cwd -P compbio_lab -l h_vmem=25G -l h_rt=5:00:00 -N binarize_p${model} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$MAINJOB --chr $chr"
done
qsub -cwd -P compbio_lab -l h_vmem=25G -l h_rt=2:00:00 -N binarize_p${model}_all -hold_jid binarize_p${model} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$MAINJOB"


# TODO: Move marks to better directories:

# Marks:
while read mark; do
    echo ${mark}
    while read chr size; do
        qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N binarize_${mark} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate pytorch_env > /dev/null; python ${BINDIR}/aggregate_binary_txt_to_npy.py --dir ${CHMM_FMTDIR}/${mark}/ --out ${CHMM_FMTDIR}/${mark}_all_bin --suffix _binary.txt.gz --mid _c_t.sub_ --chr $chr"
    done < ${CHROMSIZES_noY}

    qsub -cwd -P compbio_lab -l h_vmem=40G -l h_rt=1:30:00 -N binarize_${mark}_all -hold_jid binarize_${mark} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate pytorch_env > /dev/null; python ${BINDIR}/aggregate_binary_txt_to_npy.py --dir ${CHMM_FMTDIR}/${mark}/ --out ${CHMM_FMTDIR}/${mark}_all_bin --suffix _binary.txt.gz --mid _c_t.sub_"
done < $MARKS_LIST



# Run a series of plotting jobs:
# qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr1 -a 1000000 -b 2000000 -n 15"
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr1 -a 1000000 -b 1250000 -n 15"
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr10 -a 1000000 -b 2000000 -n 15"
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr10 -a 1000000 -b 1500000 -n 15"
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr19 -a 1000000 -b 2000000 -n 15"
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr19 -a 1000000 -b 1500000 -n 15"


qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr10 -a 0 -b 2000000 -n 15"
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr19 -a 2000000 -b 3000000 -n 15"
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr2 -a 1000000 -b 2000000 -n 15"
qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:30:00 -N plotchmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/extract_plot_chromHMM_region.sh -c chr3 -a 1000000 -b 2000000 -n 15"

