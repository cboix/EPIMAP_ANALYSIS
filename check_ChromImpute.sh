#!/bin/bash
# Check ChromImpute files:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh

IFS=$'\t'

echo "Checking CONVERTED FILES in ${CONVERTED_DIR}"
while read chr size
do
    echo $chr
    gzip -t ${CONVERTED_DIR}/${chr}_FINAL_H[1-4]*.gz >> ${CHECK_DIR}/converted_data_check_${TODAY} 2>&1
    gzip -t ${CONVERTED_DIR}/${chr}_FINAL_DNase*.gz >> ${CHECK_DIR}/converted_data_check_${TODAY} 2>&1
    gzip -t ${CONVERTED_DIR}/${chr}_FINAL_ATAC*.gz >> ${CHECK_DIR}/converted_data_check_${TODAY} 2>&1
done < $CHROMSIZES


echo "Checking FEATURE FILES in ${TRAIN_DIR}"
gzip -t ${TRAIN_DIR}/attributes_*.gz >> ${CHECK_DIR}/feature_data_check_${TODAY} 2>&1
while read chr size
do
    echo $chr
    gzip -t ${TRAIN_DIR}/${chr}_traindata_*.gz >> ${CHECK_DIR}/feature_data_check_${TODAY} 2>&1
done < $CHROMSIZES


echo "Checking PREDICTOR FILES in ${PREDICTOR_DIR}"

while read mark
do
    echo $mark
    gzip -t ${PREDICTOR_DIR}/useattributes*_${mark}_*.gz >> ${CHECK_DIR}/predictor_data_check_${TODAY} 2>&1
    gzip -t ${PREDICTOR_DIR}/classifier*_${mark}_*.gz >> ${CHECK_DIR}/predictor_data_check_${TODAY} 2>&1
done < $MARKS_LIST

# Make predictor file list (unsorted):
if [[ ! -s $CIDIR/predictor_dir_files ]];then
    cd ${PREDICTOR_DIR}
    ls -U > ${CIDIR}/predictor_dir_files
fi

# Remove all old files (currently runs in impute task - checks for older than may)
cd $PREDICTOR_DIR
for TASK in `seq 1 9919`
do
    # Get sample/mark from table (available or to_impute):
    SAMPLE=$( sed "${TASK}q;d" ${FULLIMPUTATION_TAB} | awk -v FS="\t" '{print $1}' )
    MARK=$( sed "${TASK}q;d" ${FULLIMPUTATION_TAB} | awk -v FS="\t" '{print $2}' )
    echo "TASK $TASK - checking $SAMPLE and $MARK"

    # Remove old files:
    cmd="rm $( ls -clt $( grep ${SAMPLE}_${MARK} ${CIDIR}/predictor_dir_files ) | awk '$0 !~/ May /{print $9}' | awk '{printf $1" "}' )"
    if [[ "$cmd" != "rm " ]];then
        echo "Removing files..."
        bash -c "$cmd"
    fi
done < $CIDIR/cell_types.tsv



# CHECK IMPUTED FILES AND REMOVE FILES THAT FAIL
echo "Checking IMPUTED FILES in ${IMPUTED_DIR}"
while read chr size
do
    echo $chr
    gzip -t ${IMPUTED_DIR}/${chr}_impute_*.gz >> ${CHECK_DIR}/imputed_data_check_${TODAY} 2>&1
done < $CHROMSIZES

echo "[STATUS] Removing all EOF error imputed files"
rmcmd="rm $(grep "unexpected end of file" ${CHECK_DIR}/imputed_data_check_${TODAY} | awk -vFS=": " '{print $2}' )"
if [[ "$rmcmd" != "rm " ]];then
    echo "Removing files..."
    echo $rmcmd
    bash -c "$rmcmd"
fi

