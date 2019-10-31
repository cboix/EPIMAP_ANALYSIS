#!/bin/bash 
export CHIPDIR=$DBDIR/ChIP-seq
export LNKDIR=$CHIPDIR/file_links/all_submitted_released
# ========================================
# 1. Turn each BAM file into TagAlign file
# ========================================
source $BINDIR/ChIP_download_preprocess.sh

# ==========================================
# 2. Pool data and determine fragment length
# Also determine number of reads in total: 
# Subsample if necessary! Or normalize some other way
# ==========================================
source $BINDIR/ChIP_2_Pool_and_QC.sh

# ===========================
# 3. Call peaks with controls
# 3b. Turn file into bedgraph 
# ===========================
source $BINDIR/ChIP_3_Call_Peaks_into_BG.sh

