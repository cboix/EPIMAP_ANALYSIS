#!/bin/bash
# ------------------------------------------
# Assign TSS to promoter elements,
# Make possible links using same coordinates
# ------------------------------------------
# NOTE: Run on kgpu2 for now - use for load link data

source /home/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 
# TODO: Change if want to use nonoverlapping
AGAINST=$CORECOORD

cd $DBDIR/GTEx_eQTL/GTEx_Analysis_v7_eQTL

# Get TSS:
awk -vOFS="\t" 'NR > 1{tss=($6=="-"?$5:$4); print "chr"$3, tss, tss, $1}' all_egenes.txt | sort -k1V -k2n > all_egenes_tss.tsv

# Merge with all elements (for now) - later filter for only promoters:
bedtools closest -a all_egenes_tss.tsv -b $AGAINST -d | awk -vOFS="\t" '{print $1,$2,$4,$8,$9,NR}' > tss_map_closest.tsv

# Create TSS regions + get all testable overlaps:
WINDOW=1000000  # 2MB window
awk -vOFS="\t" -v win=$WINDOW '{min=$2-win; min=(min > 0?min:0); print $1, min, $2+win, $6}' tss_map_closest.tsv > tss_1MB_regions.tsv

# Intersect to get all overlaps:
bedtools intersect -a $AGAINST -b tss_1MB_regions.tsv -wb | awk -vOFS="\t" '{print $4, $8}' | gzip -c > all_1MB_intersections.tsv.gz

# Count of overlaps:
zcat all_1MB_intersections.tsv.gz | wc -l 

# All intersections:
bedtools intersect -a all_egenes_tss.tsv -b $AGAINST -wb | cut -f4 | sort -u | wc

# bedtools intersect -a all_egenes_tss.tsv -b $AGAINST -wb | cut -f4 | sort -u | wc
cd $ANNDIR
GENCTX=$ANNDIR/gencode.tx.v30lift37.basic.annotation.gtf.gz
GENCTSS=$ANNDIR/gencode.tx.tss.v30lift37.basic.annotation.gtf.gz
zcat $GENCTX | awk -vOFS="\t" -F"\t" '$1 ~/chr[0-9X]/{gsub(";.*$","", $9); tss=($7=="-"?$5:$4); print $1, tss, tss, $9}' | sort -k1V -k2n -u | gzip -c > $GENCTSS


bedtools intersect -a $GENCTSS -b $AGAINST -wb | cut -f4 | sort -u | wc

# ----------------------
# Make validation links:
# ----------------------
# Promoter Capture HiC
# ----------------------
# TODO: GET ORIGINAL DATA - supplement Javierre
mkdir -p $DBDIR/linking_data
cd $DBDIR/linking_data

# Using Promoter Capture HiC (Javierre 2016):
PCFILE=ActivePromoterEnhancerLinks_JavierrePChIC.tsv.gz

zcat $PCFILE | awk -vOFS="\t" 'NR > 1{print $1,$2,$3,$4}' | sort -k1V -k2n -u > tmpbait.bed
zcat $PCFILE | awk -vOFS="\t" 'NR > 1{print $5,$6,$7,$8}' | sort -k1V -k2n -u > tmpoe.bed

# NOTE: bait is long regions, filter by promoter only afterwards. - captures 176k with just 8k bait tags.
bedtools intersect -a $AGAINST -b tmpbait.bed -wb | awk -vOFS="\t" 'BEGIN{print "name", "baitID"}{print $4, $8}' > baitint.tsv

# Capture also long -from 25k to 250k
bedtools intersect -a $AGAINST -b tmpoe.bed -wb | awk -vOFS="\t" 'BEGIN{print "name", "oeID"}{print $4, $8}' > oeint.tsv

# Need to decide how to filter these down?


# ----------------------
# Roadmap links:
# http://www.biolchem.ucla.edu/labs/ernst/roadmaplinking/
# ----------------------
cd $DBDIR/linking_data/RoadmapLinks

cat links_E*.txt | awk -vOFS="\t" '{print $1,$2,$3,$4}' | sort -k1V -k2n -u | gzip -c > all_links.tsv.gz

# NOTE: bait is long regions, filter by promoter only afterwards. - captures 176k with just 8k bait tags.
bedtools intersect -a $AGAINST -b all_links.tsv.gz -wb | awk -vOFS="\t" 'BEGIN{print "name", "GeneID"}{print $4, $8}' | sort -u | gzip -c > all_links_mapped.tsv.gz
# TODO: top is diff.

# About 1e6 links total:
zcat all_links_mapped.tsv.gz | wc -l 


