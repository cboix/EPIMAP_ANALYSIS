#!/bin/bash
# ------------------------------------------
# Munge the CAD example summary statistics
# ------------------------------------------
start=`date +%s`
hostname -f
source /home/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES
# source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES

cd $DBDIR/gwas_sumstats
mkdir -p munged 

# -----------------------------------------------
# 1. Re-format all of the relevant GWAS sumstats:
# -----------------------------------------------
# CAD: CAD_UKBIOBANK.gz  # hg19
zcat CAD_UKBIOBANK.gz | awk -vOFS="\t" 'NR>1{print "chr"$2,$3,$6,$8}' | sort -k1,1V -k2,2n | gzip -c > munged/CAD.munged.tsv.gz

# SBP: systolic_UKB2_sumstats.txt.gz # hg19
zcat systolic_UKB2_sumstats.txt.gz | awk -vOFS="\t" 'NR>1{print "chr"$2,$3,$7,$12}' | sort -k1,1V -k2,2n | gzip -c > munged/SBP.munged.tsv.gz

# HDL: jointGwasMc_HDL.txt.gz # Is from Willer 2013 # hg19
zcat jointGwasMc_HDL.txt.gz | awk -vOFS="\t" 'NR>1{sub(":","\t", $2); print $2,$6,$9}' | sort -k1,1V -k2,2n | gzip -c > munged/HDL.munged.tsv.gz

# AFB: nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz # hg19
zcat nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz | awk -vOFS="\t" 'NR>1{print "chr"$3,$4,$8,$10}' | sort -k1,1V -k2,2n | gzip -c > munged/AFB.munged.tsv.gz

# WAB (merge):
# WHRADJ.WOMEN.LE50.publicrelease.txt.gz
# WHRADJ.WOMEN.GT50.publicrelease.txt.gz
# WHRADJ.MEN.LE50.publicrelease.txt.gz
# WHRADJ.MEN.GT50.publicrelease.txt.gz
zcat whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz | awk -vOFS="\t" 'NR>1{sub(":.*","",$3); print $3,"chr"$1,$2}' | sort +0 -1 | gzip -c > rsmapping.tsv.gz
zcat rsmapping.tsv.gz > rsmapping.tsv
zcat WHRADJ.WOMEN.LE50.publicrelease.txt.gz | awk -vOFS="\t" 'NR>1{print $1,$5,$7}' | sort +0 -1 | join - rsmapping.tsv | awk -vOFS="\t" '{print $4,$5,$2,$3,"wl5"}' | sort -k1,1V -k2,2n | gzip -c > munged/WHR.WL5.munged.tsv.gz
zcat WHRADJ.WOMEN.GT50.publicrelease.txt.gz | awk -vOFS="\t" 'NR>1{print $1,$5,$7}' | sort +0 -1 | join - rsmapping.tsv | awk -vOFS="\t" '{print $4,$5,$2,$3,"wg5"}' | sort -k1,1V -k2,2n | gzip -c > munged/WHR.WG5.munged.tsv.gz
zcat WHRADJ.MEN.LE50.publicrelease.txt.gz | awk -vOFS="\t" 'NR>1{print $1,$5,$7}' | sort +0 -1 | join - rsmapping.tsv | awk -vOFS="\t" '{print $4,$5,$2,$3,"ml5"}' | sort -k1,1V -k2,2n | gzip -c > munged/WHR.ML5.munged.tsv.gz
zcat WHRADJ.MEN.GT50.publicrelease.txt.gz | awk -vOFS="\t" 'NR>1{print $1,$5,$7}' | sort +0 -1 | join - rsmapping.tsv | awk -vOFS="\t" '{print $4,$5,$2,$3,"mg5"}' | sort -k1,1V -k2,2n | gzip -c > munged/WHR.MG5.munged.tsv.gz

# WHR: whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz
zcat whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz | awk -vOFS="\t" 'NR>1{print "chr"$1,$2,$7,$9}' | sort -k1,1V -k2,2n | gzip -c > munged/WHR.munged.tsv.gz

# ---------------------------------------
# 2. Merge each against the CAD variants:
# ---------------------------------------
cd munged
zcat CAD.munged.tsv.gz | awk -vOFS="\t" '{print $1"_"$2,$3,$4}' | sort +0 -1 > CAD_tmp.tsv

gwabbr=( SBP HDL AFB WHR.WL5 WHR.WG5 WHR.ML5 WHR.MG5 WHR )
for gwas in ${gwabbr[@]}; do
    echo $gwas
    if [[ ! -s ${gwas}.wCAD.tsv.gz ]]; then
        zcat ${gwas}.munged.tsv.gz | awk -vOFS="\t" '{print $1"_"$2,$3,$4}' | sort +0 -1 | join - CAD_tmp.tsv | awk -vOFS="\t" '{sub("_","\t", $1); print $1,$2,$3,$4,$5}' | sort -k1,1V -k2,2n | gzip -c > ${gwas}.wCAD.tsv.gz
    fi
done


end=`date +%s`
runtime=$((end-start))
echo "Finished $COMMAND sucessfully in $runtime seconds."
