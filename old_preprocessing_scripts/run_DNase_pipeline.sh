#!/bin/bash
# SGE Array to run step 1 for a single replicate -- merge and norm afterwards.
# run=$( sed "${SID}q;d" $DHSINFO ) # FIXME only for testing
# Arguments (feed files to process:)
FILEINFO=$1
if [[ -s $FILEINFO ]]
then
    run=$( sed "${SGE_TASK_ID}q;d" $FILEINFO )
else 
    run=$( sed "${SGE_TASK_ID}q;d" $DHSINFO )
fi
eval $( echo ${run} | awk '{sub("tered align","tered_align",$0); a = $1; sub(".*download/","",a); sub("\\..*$","",a); printf("id=\"%s\"; link=\"%s\"; assembly=\"%s\"; cell=\"%s\"; processing=\"%s\";",a,$1,$3,$5,$4)}' )

echo "$run"
echo "$id $link $assembly $cell $processing"

# ========================================
# Preprocess and turn into tagAlign files:
# ========================================
echo "- STEP1 Preprocess"
STEP1_FILE="${FILE_DIR}/FINAL_${id}_${cell}.tagAlign.gz"
echo ${STEP1_FILE}
if LC_ALL=C gzip -l ${STEP1_FILE} | awk 'NR==2 {exit($2!=0)}'
then
    source $BINDIR/DNase_download_preprocess.sh $id $cell $link $assembly
fi

# ===========================================
# Create appropriate windows for the dataset:
# ===========================================
echo "- STEP2 Coverage into TSS Windows:"
export DW_FILE=${WINDIR}/${id}_${cell}.DHS_windows
echo ${DW_FILE}
if [[ ! -s $DW_FILE  || "9" != "$( awk 'NR==2{print NF}' $DW_FILE )" ]] 
then
    source $BINDIR/DNase_compute_window_occupancy.sh $id $cell
fi

# TODO May need to do these steps outside of this file:
echo "- STEP3 Pool and Subsample"
# # In serial:
# # For each epitope + DNase + WCE:
# IFS=$'\t'
# while read epitope
# do 
#     echo "- STEP2 ${cell}_${epitope}"
#     STEP2_FILE="${CELL_DIR}/FINAL_${cell}_${epitope}.tagAlign.gz"
#     echo $STEP2_FILE
#     if [[ ! -s ${STEP2_FILE} ]]
#     then
#         source $BINDIR/DNase_download_preprocess.sh $cell $epitope ${CELL_DIR}
#     fi
# done < $DBDIR/epitopes # list of epitopes we are interested in

# TODO Figure out if there are enough reads in the full pooled dataset!

echo "- STEP4 Call DHS peaks + coverage of putative Promoter and ENH elements"


# # PATCHING:
# for SID in `seq 2 8`
# do
#     run=$( sed "${SID}q;d" $DHSINFO )
#     eval $( echo ${run} | awk '{sub("tered align","tered_align",$0); a = $1; sub(".*download/","",a); sub("\\..*$","",a); printf("id=\"%s\"; link=\"%s\"; assembly=\"%s\"; cell=\"%s\"; processing=\"%s\";",a,$1,$3,$5,$4)}' )

#     echo "$run"
#     echo "$id $link $assembly $cell $processing"

#     # ===========================================
#     # Create appropriate windows for the dataset:
#     # ===========================================
#     echo "- STEP2 Coverage into TSS Windows:"
#     export DW_FILE=${WINDIR}/${id}_${cell}.DHS_windows
#     echo ${DW_FILE}
#     if [[ ! -s $DW_FILE  || "9" != "$( awk 'NR==2{print NF}' $DW_FILE )" ]] 
#     then
#         source $BINDIR/DNase_compute_window_occupancy.sh $id $cell
#     fi
# done

