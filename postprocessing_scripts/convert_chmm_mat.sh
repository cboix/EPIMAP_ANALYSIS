#!/bin/bash/
# qsub -cwd -t 1-23 -l h_vmem=20G -l h_rt=4:00:00 -N split_chmm -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source $BINDIR/convert_chmm_mat.sh"

source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
chr=$( cut -f1 $CHROMSIZES_noY | sort -u | sed "${SGE_TASK_ID}q;d" )
echo "$chr"

CDATADIR=/broad/compbio_ce/cboix/EPIMAP_ANALYSIS/updated_hg38_chromHMM
mkdir -p $CDATADIR/ $CDATADIR/split

cd $CDATADIR

echo "[STATUS] Running split chrom:"
t1=`date +%s`

NAMFILE=$DBDIR/public_metadata_released/samplelist_833_final.tsv
zcat calls_${chr}_hg38.bed.gz | awk -v chr=$chr -vOFS="\t" -F"\t" 'BEGIN{a=0} FNR==NR {a=a+1; bss[a]=$1; next} {for(i=4;i<=NF;i++){print $1,$2,$3,"E"$i > "split/"bss[i-3]"_"chr"_segments.bed"}}' $NAMFILE -

t2=`date +%s`
runtime=$((t2-t1))
echo "Finished $chr in $runtime seconds."

