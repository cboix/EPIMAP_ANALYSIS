#!/bin/bash
#SBATCH -J process_IC_files        # Job name (-N)
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G              # Job memory request (-l h_vmem)
#SBATCH --time=01:00:00         # Time limit hrs:min:sec (-l h_rt)
#SBATCH -p kellis
#SBATCH --output=/home/cboix/data/EPIMAP_ANALYSIS/db/out/IC/process_IC_files_%j.log
#SBATCH --export=ALL            # Copy env. variables (-V)
#SBATCH --requeue               # For if requeue needed (-r y)

# Grid Engine options
#$ -N process_IC_files
#$ -cwd
# #$ -P compbio_lab
#$ -l h_vmem=10G 
#$ -l h_rt=01:00:00
#$ -tc 250
#$ -j y
#$ -b y 
#$ -V 
#$ -r y 
#$ -o /broad/compbio/cboix/EPIMAP_ANALYSIS/db/out/IC
#$ -e /broad/compbio/cboix/EPIMAP_ANALYSIS/db/out/IC
#$ -t 1-1428

# 1428 runs:
# sbatch --array=1-100 $BINDIR/ic_eval/process_IC_files.sh
# qsub -t 1-100 $BINDIR/ic_eval/process_IC_files.sh

# Config and directories:
# -----------------------
source ${HOME}/data/EPIMAP_ANALYSIS/bin/ic_eval/config_IC_files.sh

# Cluster type: 
if [[ "$SLURM_JOB_ID" == "" ]]; then
    TASK=$SGE_TASK_ID 
else 
    TASK=$SLURM_ARRAY_TASK_ID
fi

# ------------------------------------------------
# Figure out which track to process for each task:
# ------------------------------------------------
nteam=$( wc -l $TEAMFILE | cut -f1 -d" " )
ntrack=$( wc -l $INFOFILE | cut -f1 -d" " )
let "ngrp=1+$nteam"

# Track:
let "tracknum=$TASK % $ntrack"
let "grpnum=($TASK - $tracknum)/ $ntrack"
sample=$( awk -v id=$tracknum 'NR==(id+1){print $1}' ${INFOFILE} )
mark=$( awk -v id=$tracknum 'NR==(id+1){print $2}' ${INFOFILE} )
track=$( awk -v id=$tracknum 'NR==(id+1){print $3}' ${INFOFILE} )

# Group (or blind):
if [[ "$grpnum" -eq 0 ]]; then
    num="blind"
    team="blind"
elif [[ "$grpnum" -eq 27 ]]; then
    num="qnorm"
    team="qnorm"
else
    num=$( awk -v id=$grpnum 'NR==id{print $1}' ${TEAMFILE} )
    team=$( awk -v id=$grpnum 'NR==id{print $2}' ${TEAMFILE} )
fi

# -------------------------------
# Download and process the track:
# -------------------------------
export TMP_DIR=$TMP/${sample}${mark}_process_${RANDOM}
mkdir -p ${TMP_DIR}
echo "[STATUS] Downloading + processing: $sample $mark $track for team $team"

cd $ICDIR

# Download directory and parent URL where file is:
# ------------------------------------------------
if [[ "$team" == "blind" ]]; then
    parentdir=${baseurl}/blind
    DWDIR=${BLDIR}
elif [[ "$team" == "qnorm" ]]; then
    parentdir=${qnormurl}
    track=${sample}${mark}.qnorm.bigwig
    DWDIR=${BLDIR}
else
    if [[ $num -lt 101 ]]; then
        # Avocado + Average
        lteam=$( echo "$team" | tr '[:upper:]' '[:lower:]' )
        parentdir=${baseurl}/${lteam}
        track=${sample}${mark}.bigwig
    elif [[ $num -eq 200 ]]; then
        # ChromImpute from Jason
        parentdir=${ciurl}
        track=${sample}_${mark}.bw
    else
        # Team predictions:
        parentdir=${baseurl}/round2/${num}
        track=${sample}${mark}.bigwig
    fi
    DWDIR=${R2DIR}/${num}
fi

# Files to create:
# ----------------
export BWFILE=${DWDIR}/${num}_${sample}${mark}.bigwig
export BDGFILE=${DWDIR}/${num}_${sample}${mark}.bedgraph
export WIGFILE=$WIGDIR/chr21_${num}_${sample}${mark}.bedgraph.wig.gz
mkdir -p ${DWDIR}

if [[ -s $WIGFILE ]]; then 
    out=$( gzip -t $WIGDIR/*_${num}_${sample}${mark}.bedgraph.wig.gz ) 
    if [[ "$out" != "" ]]; then
        echo "Removing all wig files"
        rm $WIGDIR/*_${num}_${sample}${mark}.bedgraph.wig.gz
    fi
fi

# Download the bigwig file:
# -------------------------
if [[ ! -s $WIGFILE ]]; then
    if [[ ! -s $BWFILE ]] && [[ ! -s $BDGFILE ]]; then
        echo "[STATUS] Downloading track: ${parentdir}/${track}"
        if [[ $num -eq 200 ]]; then
            wget ${parentdir}/${track} -O $BWFILE --no-check-certificate
        else
            wget ${parentdir}/${track} -O $BWFILE
        fi
    else
        echo "[STATUS] Bigwig ${track} already exists for team ${team} (${num})"
        ls -sh $BWFILE
    fi

    # Convert track to bedgraph:
    # --------------------------
    if [[ ! -s $BDGFILE ]]; then
        echo "[STATUS] Converting track to bedgraph: ${track} for team ${team}"
        bigWigToBedGraph $BWFILE $BDGFILE
    else
        echo "[STATUS] Bedgraph ${track} already exists for team ${team} (${num})"
        ls -sh $BDGFILE
    fi

    # Convert track to bedgraph:
    # --------------------------
    if [[ ! -s $WIGFILE ]]; then
        echo "[STATUS] Converting track to wig files for Eval: ${track} for team ${team}"
        export CVFILE=${TMP_DIR}/convert_info.tsv
        echo -e "$sample\t$mark\t${num}_${sample}${mark}.bedgraph" > $CVFILE
        # ls -sh $CVFILE
        cat $CVFILE
        cmd="$CIBASE Convert ${DWDIR} $CVFILE $CHROMFILE ${WIGDIR}"
        echo "$cmd"
        bash -c "${cmd}"
        # java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar $CHROMIMPUTE Convert ${DWDIR}/ $CVFILE $CHROMFILE ${WIGDIR}/
    else
        echo "[STATUS] Wig files for ${track} already exist for team ${team} (${num})"
        ls -sh $WIGDIR/chr*_${num}_${sample}${mark}.bedgraph.wig.gz
    fi

    # Clean up files
    if [[ -s $WIGFILE ]]; then
        rm $BWFILE $BDGFILE
    fi
else
    echo "[STATUS] Wig files for ${track} already exist for team ${team} (${num})"
    ls -sh $WIGDIR/chr*_${num}_${sample}${mark}.bedgraph.wig.gz
fi

rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished track download and conversion sucessfully in $runtime seconds."
