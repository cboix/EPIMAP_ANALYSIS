#!/bin/bash
# Process GTEx eQTLs - python overlap modules are not functional.

source /home/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES
# source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES

# export TMP_DIR=${TMP}/CTOM_${COMMAND}_${MODEL}_${RANDOM}
# mkdir -p ${TMP_DIR}

cd $DBDIR/GTEx_eQTL/GTEx_Analysis_v7_eQTL

# 1. Create SNP list:
rm tmpfile
while read file; do
    echo $file
    zcat $file | awk 'NR > 1{split($1,a, "_"); print "chr"a[1],a[2]}' >> tmpfile
done < <( ls *.v7.signif_variant_gene_pairs.txt.gz)


# 2. Reduce file - from 36M to 3M
sort -k1V -k2n -u tmpfile | awk -vOFS="\t" '{print $1, $2, $2}' > snplist
wc -l tmpfile
wc -l snplist

head -n 5 snplist 
rm tmpfile

# 3. Get all overlapping or closest chunks:
bedtools closest -a snplist -b $CORECOORD -d | awk -vOFS="\t" '{print $1,$2,$7,$8}' > int_snplist.tsv

# All snps
cut -f1,2 int_snplist.tsv | sort -u | wc -l 

# Reduced sets close to enhancers:
awk '$4 <= 500{print $1,$2}' int_snplist.tsv | sort -u | wc -l 
awk '$4 <= 100{print $1,$2}' int_snplist.tsv | sort -u | wc -l 
awk '$4 <= 0{print $1,$2}' int_snplist.tsv | sort -u | wc -l 
