#!/bin/bash
# --------------------------------------
# Merge peaks into binary matrix, 
# 3. Create a final 0/1 file as:
#    a. Sparse in MatrixMarket fmt (mtx)
#    b. Pickled CSR Matrix for python
# --------------------------------------
MARK=$1
PREF=$2
TMP_DIR=$3

# Prefixes:
OUTPREFIX=${MPD_DIR}/master_${MARK}_from${PREF}
MTXTMP=${TMP_DIR}/tmp_assign.mtx
# 1. Concatenate all DHS:
cat ${TMP_DIR}/DHSs_all/chunk*.bed > ${OUTPREFIX}_DHSs_info.tsv

# 2. Concatenate all peak observations:
rm $MTXTMP
while read file; do
    echo $file
    awk -vOFS="\t" -v pref="FINAL_${MARK}_" '{sub("impute_","imputed\t",$4); sub(pref,"observed\t",$4); sub("_.*","",$4); sub(".sub","",$4); print $8,$4,$5}' $file >> $MTXTMP
done < <( ls ${TMP_DIR}/peaks_all/chunk*bed | sort -V )

# Assign numbers for each dhs loc (mapping pks)
cut -f1 $MTXTMP | uniq | awk -vOFS="\t" '{print $1,NR}' > ${TMP_DIR}/loc_list.tsv

# Assign num for each possible CT (mapping ct)
cut -f3 $MTXTMP | uniq | sort -u | awk -vOFS="\t" '{print $1,NR}' > ${TMP_DIR}/cell_list.tsv

# Size of matrix:
ROWNUM=$( cat $TMP_DIR/loc_list.tsv | wc -l )
COLNUM=$( cat $TMP_DIR/cell_list.tsv | wc -l )
echo "[STATUS] Have $ROWNUM locations in $COLNUM cells."

echo "[STATUS] Processing tmp matrix into MTX and CSR formats"
python ${BINDIR}/proc_mtx_to_sparse.py --mtx $MTXTMP --cells ${TMP_DIR}/cell_list.tsv --out $OUTPREFIX  # --nrow $ROWNUM --ncol $COLNUM

# Look at outfiles:
ls -sh ${OUTPREFIX}*

if [[ -s ${MPD_DIR}/${OUTPREFIX}_imputed_csr.cp ]]; then
    echo "[STATUS] CSR matrices created. Removing temp dir."
    rm -rf ${TMP_DIR}
fi

end=`date +%s`
runtime=$((end-start))
echo "Finished submitting masterlist for $MARK with prefix: $PREF using file: $INFOFILE in $runtime seconds."
