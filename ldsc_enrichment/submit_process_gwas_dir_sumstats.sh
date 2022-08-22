#!/bin/bash
# -----------------------------------
# Turn all tabix files in GWAS dir to 
# sumstats.gz format for LDSC
# -----------------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

# Directories/Files:
export LDSCDIR=$SFTDIR/ldsc
export MAINGWAS=/broad/compbio/data/gwas
export H3SNPS=/broad/compbio/cboix/genomes/hapmap3/w_hm3.snplist
export GWASLIST=${DBDIR}/gwas_list.tsv
ls $MAINGWAS/*/tabix/*.bed.gz | awk -vOFS="\t" '{a=$1; split($1,arr,"/"); sub(".bed.gz","",arr[8]); print $1, arr[6], arr[8]}' > $GWASLIST

GNUM=$( cat $GWASLIST | wc -l )

echo "[STATUS] Submitting process $GNUM GWAS files to sumstats"
qsub -cwd -t 1-$GNUM -P compbio_lab -l h_vmem=5G -l h_rt=0:30:00 -N proc_gwas_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "${BINDIR}/ldsc_enrichment/process_gwas_dir_sumstats.sh $GWASLIST"


# Also create pruned leadsnp files:
qsub -cwd -t 1-$GNUM -P compbio_lab -l h_vmem=5G -l h_rt=0:15:00 -N proc_gwas_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "${BINDIR}/ldsc_enrichment/process_gwas_dir_leadsnps.sh $GWASLIST"

# Aggregate the lead snps:
GWLDIR=$DBDIR/gwas_leadsnps
mkdir -p $GWLDIR

LEADFILE=$GWLDIR/all.lead.tsv
if [[ ! -s ${LEADFILE}.gz ]]; then
    rm $LEADFILE
    while read ff subdir trait; do 
        sumpref=${GWLDIR}/${subdir}/${trait}
        file=${sumpref}.lead.tsv.gz
        echo $(basename $file)
        zcat $file | awk -vOFS="\t" -v sset=$subdir -v tt=$trait '{print $1,$2,$3,sset,tt}' >> $LEADFILE
    done < $GWASLIST
    gzip -f $LEADFILE
fi

