#!/bin/bash
# --------------------------------------------------------------------------
# Convert the Imputation Challenge data to hg19, 25bp fixed step, and prune:
# --------------------------------------------------------------------------
# conda activate mv_env + run as: 
# qsub -cwd -t 1-51 -l h_vmem=20G -l h_rt=4:00:00 -N process_ev -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/convert_IC_data_hg19_fixedstep.sh"
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
conda activate mv_env

ICDIR=/broad/compbio_ce/cboix/imputation_challenge
DATADIR=${ICDIR}/blind_data
INFOFILE=${DATADIR}/ic_prefixes.txt
CVDIR=${DATADIR}/converted
EVDIR=${DBDIR}/external_validation
mkdir -p $DATADIR ${CVDIR} $EVDIR

# Add the libpng12 library:
LD_LIBRARY_PATH="/broad/compbio/cboix/software/miniconda3/envs/mv_env/lib/:$LD_LIBRARY_PATH"
LD_LIBRARY_PATH="/broad/compbio/cboix/software/MATLAB/MATLAB_Compiler_Runtime/v714/bin/glnxa64:$LD_LIBRARY_PATH"

echo "$SGE_TASK_ID"
prefix=$( cut -f1 ${INFOFILE} | sed "${SGE_TASK_ID}q;d" )
id=$( cut -f2 ${INFOFILE} | sed "${SGE_TASK_ID}q;d" )
mark=$( cut -f3 ${INFOFILE} | sed "${SGE_TASK_ID}q;d" )
SAMPLE=$( grep $id ${DATADIR}/ic_namespace_mapping.tsv | awk '{print $1}' )
MARK=$( grep $mark ${DATADIR}/ic_mark_mapping.tsv | awk '{print $2}' )

echo "[STATUS] ${prefix} for ${id} in ${mark}"
echo "[STATUS] ${prefix} matches to ${SAMPLE} in ${MARK}"

# -------------------------
# Convert to hg38 bedgraph:
# -------------------------
bwfile=${DATADIR}/${prefix}.pooled.bigwig
bdgfile=${DATADIR}/${prefix}.pooled.bedgraph
if [[ ! -s ${bdgfile}.gz ]]; then
    if [[ ! -s $bdgfile ]]; then
        bigWigToBedGraph ${bwfile} ${bdgfile}
    fi
    gzip ${bdgfile}
fi

# ---------------------------------------------
# Convert to fixed step 200bp or 25bp bedgraph:
# ---------------------------------------------
step=200
steplist=(25 200)
for step in ${steplist[@]}; do
    echo $step
    chrom=chr1
    BCFILE=${DATADIR}/${prefix}.${chrom}.b${step}.bed.gz
    if [[ ! -s ${BCFILE} ]] || [[ $bdgfile.gz -nt $BCFILE ]]; then
        echo "[STATUS] Converting bedgraph to fixed-step"
        zcat ${bdgfile}.gz | awk -f $BINDIR/convert_to_fixed_step_bed.awk -vOFS="\t" -v step=$step -v chrom=$chrom | gzip -c > ${BCFILE}
    fi

    # -----------------------
    # Liftover the chr1 file:
    # -----------------------
    OUTBED=${DATADIR}/${prefix}.chr1.b${step}.lo.bed
    if [[ ! -s ${OUTBED}.gz ]] || [[ $BCFILE -nt ${OUTBED}.gz ]]; then
        TMPMAP=${prefix}.chr1.b${step}.unmapped.bed
        echo "[STATUS] Lifting over $BCFILE to $OUTBED"
        liftOver -minMatch=0.9 ${BCFILE} ${LOCHAIN} ${OUTBED} ${TMPMAP}
        gzip -f $OUTBED
    fi

    # -----------------------------------------
    # Fill and re-bin the liftedOver chr1 file:
    # -----------------------------------------
    FINALBED=${DATADIR}/${prefix}.chr1.b${step}.lo.final.bed.gz
    if [[ ! -s ${FINALBED} ]] || [[ $OUTBED.gz -nt $FINALBED ]]; then
        echo "[STATUS] Fixing and re-binning fixed-step bedgraph"
        zcat ${OUTBED}.gz | awk -v chrom=$chrom '$1 == chrom' | sort -k1,1V -k2,2n | bedtools merge -i - -c 4 -o mean | awk -vOFS="\t" 'BEGIN{s=0}{
                start=($2 - 1);
                if (s != start){print $1, s, start, 0};
                print $1, $2-1, $3,$4;
                s=$3;
            }' | awk -f $BINDIR/convert_to_fixed_step_bed.awk -vOFS="\t" -v step=$step -v chrom=$chrom | gzip -c > $FINALBED
    fi

    # -----------------------------------------------------------------------------
    # Check for gaps in the final file:
    # NOTE: Doesn't end completely, but thats ok. shouldn't have any internal gaps.
    # -----------------------------------------------------------------------------
    GAPS=$(zcat ${FINALBED} | awk 'BEGIN{end=0}{if (end != $2 - 1){print end, $0}; end=$3}' | wc -l) 
    if [ $GAPS -gt 0 ]; then
        echo "[WARNING] Greater than 0 gaps: ${GAPS} found"
    else
        echo "[STATUS] No gaps found, copying to validation analysis directory."
        cp ${FINALBED} ${EVDIR}
    fi
done


# --------------------------------------------------------
# Pull up the matching binarized and raw imputed datasets:
# --------------------------------------------------------
CBINDIR=$CIDIR/binarized_imputed/${MARK}
IMPBIN=${CBINDIR}/${SAMPLE}_${chrom}_binary.txt.gz
IMPFILE=${IMPUTED_DIR}/${chrom}_impute_${SAMPLE}_${MARK}.wig.gz
ls -sh $IMPBIN
ls -sh $IMPFILE

zcat ${IMPBIN} | awk 'NR > 2' | gzip -c > ${EVDIR}/${prefix}.${chrom}.imputed.binary.txt.gz
zcat ${IMPFILE} | awk -v count=8 'BEGIN{a=0;i=0;}NR > 2{i=i+1; a=a+$1; if(i==count){print a / count; a=0; i=0}}' | gzip -c > ${EVDIR}/${prefix}.${chrom}.imputed.b200.raw.txt.gz
zcat ${IMPFILE} | awk 'NR > 2' | gzip -c > ${EVDIR}/${prefix}.${chrom}.imputed.b25.raw.txt.gz


# ----------------------
# Calculate the metrics:
# ----------------------
# NOTE: We won't evaluate locations that are exactly 0, these are mis-mapped loc
TMPINFO=${TMP}/tmpinfo_${mark}_${MARK}_${RANDOM}.tsv
awk -v mark=$MARK '$2 == mark && $3 ~ /FINAL/' ${ALL_TRACKS_TAB} > ${TMPINFO}

if [[ ! -s ${EVDIR}/${prefix}.chr1.stats.tsv ]]; then
    $SFTDIR/miniconda3/envs/mv_env/bin/R --slave -f $BINDIR/compute_IC_data_metrics.R --args ${EVDIR}/${prefix} ${TMPINFO}
fi

