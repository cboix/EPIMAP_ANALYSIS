#!/bin/bash
# -----------------------------------
# Establish the following analysis pipeline:
# 1. Clustering of promoters and enhancers
# - Assessing methods for clustering (to do after pipeline is complete)
#
# 2. GREAT enrichment of clusters 
# 3. Motif enrichment of clusters
# 4. GWAS enrichment of clusters
# 5. GWAS enrichment of H3K27ac peaks
# 6. GWAS enrichment on a tree (with Yongjin)
# 
# 7. Linking of promoters and enhancers through their clusters and distance
# - Evaluation of linking against eQTLs, HiC
#
# TODO: In another:
# 8. Average epigenomes + their enrichments:
# - For each tissue group
# - For all selected nodes (see tree figure)
# -----------------------------------
NUMSTATES=18
NCLUST=300
ELEM="ENH"
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh $NUMSTATES

# 1. Clustering prom/enh: (preprocess, aggregate, matrix, cluster)
# TODO: SAVE ALL FINAL MATRICES TO UNIFORM DIRECTORY
# TODO: SAVE ALL FINAL CLUSTERINGS TO UNIFORM DIRECTORY
$BINDIR/submit_chmm_to_modules.sh


# 2. GREAT enrichment of clusters 
$BINDIR/submit_enrichment_pipeline.sh ENH GREAT


# 3. Motif enrichment of clusters
$BINDIR/submit_enrichment_pipeline.sh ENH MOTIF


# 4. GWAS enrichment of clusters
$BINDIR/submit_enrichment_pipeline.sh ENH GWAS


# 5. GWAS enrichment of H3K27ac peaks
# TODO: Of matrix or of peak calls?
$BINDIR/submit_enrichment_pipeline.sh H3K27ac GWAS


# 6. GWAS enrichment on a tree (with Yongjin)


# 7. Linking of promoters and enhancers through their clusters and distance
# - Evaluation of linking against eQTLs, HiC


# TODO: In another:
# 8. Average epigenomes + their enrichments:
# - For each tissue group

