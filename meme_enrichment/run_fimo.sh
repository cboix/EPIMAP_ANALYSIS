#!/bin/bash
# ----------
# FIMO Pipe:
# ----------
# 0. Load in variables:
COORDFILE=$CORECOORD
DBRES=${DML_DIR}/fimo_counts_redundant
TASK=${SGE_TASK_ID}
USEDB=0

if [ $# -eq 0 ]; then
    echo >&2 "USAGE: $(basename $0) -i [INFOFILE] -o [RESULTDIR] (OTHERARGS)
    -i  Infofile [REQUIRED] 
            Must include BED Filename and FILEID for the dataset
    -o  Output directory for results [REQUIRED]
    -d  Use pre-processed db (directory as given) [default: $USEDB]
    -c  DB coordinates file (for above -d command) [default: $COORDFILE]
    -r  DB results (for above -d command) [default: $DBRES]
    -t  Task (OPTIONAL)
            line from file to run. Otherwise runs SGE_TASK_ID: $SGE_TASK_ID"
    exit 1
fi

while getopts i:o:d:c:r:t: o
do      case "$o" in
    i)		INFOFILE="$OPTARG";;
    o)      RESULTDIR="$OPTARG";;
    t)		TASK="$OPTARG";;
    d)      USEDB="$OPTARG";;
    c)      COORDFILE="$OPTARG";;
    r)      DBRES="$OPTARG";;
    [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
esac
done

echo "[STATUS] INFOFILE: $INFOFILE"
echo "[STATUS] RESULTDIR: $RESULTDIR"
echo "[STATUS] TASK: $TASK"

# NOTE: Load after case catches variables:
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18
start=`date +%s` # For timing the pipeline

# Directories 
TMP_DIR=${TMP}/run_FIMO_${TASK}_${RANDOM}
mkdir -p ${TMP_DIR} ${RESULTDIR}

# Files and attributes:
INFILE=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $1}' )
FILEID=$( sed "${TASK}q;d" ${INFOFILE} | awk -v FS="\t" '{print $2}' )

# Temporary files:
TMP_FILE="${TMP_DIR}/regions_${RANDOM}.bed"
TMP_FASTA="${TMP_DIR}/regions_${RANDOM}.fa"
TMP_COUNTS="${TMP_DIR}/regions_${FILEID}_fimo.counts.tsv"
TMP_Q10="${TMP_DIR}/regions_${FILEID}_fimo_q10.counts.tsv"
FINAL_ALL="${RESULTDIR}/regions_${FILEID}_fimo.counts.tsv.gz"
FINAL_Q10="${RESULTDIR}/regions_${FILEID}_fimo_q10.counts.tsv.gz"

if [[ ! -s $FINAL_ALL ]]; then
    if [[ "$USEDB" != "0" ]]; then
        echo "[STATUS] Intersecting with pre-processed regions in db: $DBRES"
        # Get the row numbers of regions:
        awk '{print $1"_"$2"_"$3}' $INFILE | sort +0 -1 > ${TMP_DIR}/tmpcoll.tsv
        awk '{print $1"_"$2"_"$3, NR}' $COORDFILE | sort +0 -1 > ${TMP_DIR}/tmpcoord.tsv
        join ${TMP_DIR}/tmpcoll.tsv ${TMP_DIR}/tmpcoord.tsv | awk '{print $2}' | sort +0 -1 > $TMP_DIR/tmp_merge.tsv

        # Merge with pre-computed counts:
        zcat $DBRES/all_fimo.counts.lex.tsv.gz | join - ${TMP_DIR}/tmp_merge.tsv > ${TMP_DIR}/all_motifs.tsv
        zcat $DBRES/all_fimo_q10.counts.lex.tsv.gz | join - ${TMP_DIR}/tmp_merge.tsv > ${TMP_DIR}/q10_motifs.tsv

        # Aggregate the merges for final counts 
        NREG=$( cat $TMP_DIR/tmp_merge.tsv | wc -l )
        MOTIFMAP=${DBRES}/motif_id_mapping.tsv
        awk '{print $2}' ${TMP_DIR}/all_motifs.tsv | sort +0 -1 | uniq -c | awk '{print $2, $1}' > $TMP_DIR/all_collapsed.tsv
        awk '{print $2}' ${TMP_DIR}/q10_motifs.tsv | sort +0 -1 | uniq -c | awk '{print $2, $1}' > $TMP_DIR/q10_collapsed.tsv

        # Add number of locations + the motif name and database:
        awk '{print $2,$1}' $MOTIFMAP | sort +0 -1 | join - ${TMP_DIR}/all_collapsed.tsv | awk -vOFS="\t" -v uid=$FILEID -v nreg=$NREG '{md=$2; sub(".*_","",md); print $2, $3, nreg, uid, md}' | gzip -c > $FINAL_ALL
        awk '{print $2,$1}' $MOTIFMAP | sort +0 -1 | join - ${TMP_DIR}/q10_collapsed.tsv | awk -vOFS="\t" -v uid=$FILEID -v nreg=$NREG '{md=$2; sub(".*_","",md); print $2, $3, nreg, uid, md}' | gzip -c > $FINAL_Q10

    else
        # NOTE: Label each peak by peak NR and uid:
        # (Otherwise, get huge results)
        if [[ $INFILE =~ \.gz$ ]]; then
            zcat ${INFILE} | awk -vOFS="\t" -v uid=$FILEID '$1 ~ /chr[0-9X]/{sub("chr","",$1); print $1, $2, $3, uid"_"NR}' > ${TMP_FILE}
        else
            awk -vOFS="\t" -v uid=$FILEID '$1 ~ /chr[0-9X]/{sub("chr","",$1); print $1, $2, $3, uid"_"NR}' ${INFILE} > ${TMP_FILE}
        fi

        # TODO: Pad the ranges?
        # Turn data to fasta:
        echo "[STATUS] Getting sequences from file:"
        fastacmd="bedtools getfasta -fi ${REF_GENOME_FULL} -bed ${TMP_FILE} > ${TMP_FASTA}"
        echo "[STATUS] (Running) $fastacmd"
        bash -c "$fastacmd"
        ls -sh ${TMP_FASTA}

        # Locations of MEME tools/data
        MOTIFDB_DIR=$GENOMEDIR/meme_db/motif_databases/
        FIMOPATH=$SFTDIR/meme/bin/fimo
        # Motif databases to use:
        MDBLIST=(JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme 
            CIS-BP/Homo_sapiens.meme
            HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
            EUKARYOTE/jolma2013.meme)

        # Dirs for tmp files:
        mkdir -p ${RESULTDIR}/JASPAR/ $RESULTDIR/CIS-BP/ $RESULTDIR/HUMAN/ $RESULTDIR/EUKARYOTE/
        # output # of regions total
        NREGIONS=$( wc -l $TMP_FILE | awk '{print $1}' )
        echo "nregions: $NREGIONS"
        # Run fimo:
        rm $TMP_COUNTS
        rm $TMP_Q10
        for MD in ${MDBLIST[@]}; do 
            echo "[STATUS] Intersecting with DB: $MD"
            # For all regions - collate clusters.
            MDCOUNT=$RESULTDIR/${MD}_id_${FILEID}.tsv
            if  [[ ! -s $MDCOUNT ]]; then
                echo "[STATUS] Computing $MDCOUNT"
                ${FIMOPATH} --verbosity 1 -oc ${TMP_DIR} $MOTIFDB_DIR/$MD ${TMP_FASTA}
                mv ${TMP_DIR}/fimo.tsv $MDCOUNT
            fi
            # Concatenate all:
            awk -F"\t" 'NR > 1 && $1 !~ /#/ && $1 != ""{print $1}' ${MDCOUNT} | sort | uniq -c | awk -vOFS="\t" -v nr=$NREGIONS -v uid=$FILEID -v md=$MD '{print $2,$1,nr,uid,md}' > ${TMP_DIR}/fimo_counts.tsv
            # Aggregate counts past certain cutoff
            awk -F"\t" 'NR > 1 && $1 !~ /#/ && $1 != "" && $9 < 0.1{print $1}' ${MDCOUNT} | sort | uniq -c | awk -vOFS="\t" -v nr=$NREGIONS -v uid=$FILEID -v md=$MD '{print $2,$1,nr,uid,md}' > ${TMP_DIR}/fimo_q10.tsv
            cat ${TMP_DIR}/fimo_counts.tsv >> $TMP_COUNTS
            cat ${TMP_DIR}/fimo_q10.tsv >> $TMP_Q10
        done

        # Put in motif dir
        mv ${TMP_COUNTS} $RESULTDIR
        mv ${TMP_Q10} $RESULTDIR
        if [[ -s $FINAL_ALL ]]; then
            echo "[STATUS] Removing following intermediates:"
            ls ${RESULTDIR}/*/*meme_id_${FILEID}.tsv
            rm ${RESULTDIR}/*/*meme_id_${FILEID}.tsv
        fi
    fi
fi
    # Clear dir:
rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished FIMO motif counting in $runtime seconds."
