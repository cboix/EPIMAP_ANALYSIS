#!/bin/bash
# ------------------------------------
# Submit RNA-seq pipeline:
# 1. Get RNA-seq from ENCODE:
# 2. Aggregate RNA-seq to table
# 3. Normalize RNA-seq and make tables
# ------------------------------------
start=`date +%s`
hostname -f

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

# Directories:
export RNADIR=$DBDIR/RNA-seq
export RIDIR=$RNADIR/file_links/all_submitted_released/
export RNAINFO=$RIDIR/RNA-seq_tsv.csv
export FDIR=$RNADIR/files/RNA-seq
mkdir -p $RNADIR/files $FDIR $FDIR/tsv $FDIR/qc

RNUM=$(cat $RNAINFO | wc -l)

# ---------------------------
# 1. Get RNA-seq from ENCODE:
# ---------------------------
qsub -cwd -t 2-$RNUM -P compbio_lab -tc 2500 -l h_vmem=2G -l h_rt=0:15:00 -N dwrna -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_get_proc_RNAseq.sh"


# -----------------------------
# 2. Aggregate RNA-seq to table 
# -----------------------------
# TODO: Add hold from before:
echo "[STATUS] Aggregating $(ls $FDIR/tsv/*_fpkm.tsv | wc -l) files"
# Different versions have different genes:
wc $FDIR/tsv/*_fpkm.tsv | sort -n | awk 'BEGIN{a=0}$4 ~ /tsv/{if($1!=a){a=$1; print $1, $4}}' > $FDIR/tsv/rep_files.tsv

# Get the common denominator of genes:
colfile=${FDIR}/tsv/gene_columns.tsv
rm $colfile
while read ngenes file; do 
    echo $ngenes $file
    awk '{sub("\\.[0-9]*", "", $1); print $1}' $file > ${colfile}.tmp
    if [[ ! -s $colfile ]]; then
        mv ${colfile}.tmp $colfile
    else
        join ${colfile} ${colfile}.tmp > ${colfile}.tmp2
        mv ${colfile}.tmp2 $colfile
        rm ${colfile}.tmp
    fi
done < $FDIR/tsv/rep_files.tsv
echo "[STATUS] Will aggregate to $( cat $colfile | wc -l ) genes"

# Aggregate:
allfile=${FDIR}/tsv/all_expression.tsv
rm $allfile
while read file; do
    base=$(basename $file)
    ftag=${base%%_fpkm.tsv}
    echo $ftag
    # Merge + join:
    awk '{sub("\\.[0-9]*","",$1); print $0}' $file | join $colfile - | awk -v ftag=$ftag -vOFS="\t" '{print $2,ftag}' >> ${allfile}
done < <( ls $FDIR/tsv/*_fpkm.tsv )
gzip -f ${allfile}

# TODO: Aggregate differently if gene vs. transcript?
# NOTE: Move to kgpu for this next analysis (for now)


# ------------------------------------
# 3. Normalize RNA-seq and make tables
# ------------------------------------


end=`date +%s`
runtime=$((end-start))
echo "Finished RNA-seq download, process, aggregate, and normalize sucessfully in $runtime seconds."
