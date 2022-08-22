#!/bin/bash
#SBATCH -J score_IC_tracks        # Job name (-N)
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=25G              # Job memory request (-l h_vmem)
#SBATCH --time=04:00:00         # Time limit hrs:min:sec (-l h_rt)
#SBATCH -p kellis
#SBATCH --output=/home/cboix/data/EPIMAP_ANALYSIS/db/out/IC/score_IC_tracks_%j.log
#SBATCH --export=ALL            # Copy env. variables (-V)
#SBATCH --requeue               # For if requeue needed (-r y)

# Grid Engine options
#$ -N score_IC_tracks
#$ -cwd
# #$ -P compbio_lab
#$ -l h_vmem=25G 
#$ -l h_rt=04:00:00
#$ -tc 250
#$ -M cboix@mit.edu 
#$ -m a 
#$ -j y
#$ -b y 
#$ -V 
#$ -r y 
#$ -o /broad/compbio/cboix/EPIMAP_ANALYSIS/db/out/IC
#$ -e /broad/compbio/cboix/EPIMAP_ANALYSIS/db/out/IC
#$ -t 1-1326


# Config and directories:
# -----------------------
source ${HOME}/data/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

# SLURM_JOB_ID='somestring'
# SLURM_ARRAY_TASK_ID=1

# Cluster type: 
if [[ "$SLURM_JOB_ID" == "" ]]; then
    TASK=$SGE_TASK_ID 
else 
    TASK=$SLURM_ARRAY_TASK_ID
fi

export ICDIR=${DBDIR}/imputation_challenge
export BLDIR=${ICDIR}/blind
export R2DIR=${ICDIR}/round2
export STATDIR=${ICDIR}/stats
export WIGDIR=${ICDIR}/wig
mkdir -p ${ICDIR} ${BLDIR} ${R2DIR} ${STATDIR} ${WIGDIR}

# Files for run:
export baseurl="http://mitra.stanford.edu/kundaje/ic"
export ciurl="https://public.hoffman2.idre.ucla.edu/ernst/C86ZP"
export CIBASE="java -mx800M -jar $BINDIR/ChromImpute_SRC/ChromImpute_coeffv.jar"
export INFOFILE=${ICDIR}/ic_blind_tracks.tsv
export TEAMFILE=${ICDIR}/team_name_round2.tsv
if [[ ! -s $TEAMFILE ]]; then
    wget ${baseurl}/round2/team_name_round2.tsv -O $TEAMFILE
    echo "200\tChromImpute" >> $TEAMFILE
fi

# Eval files:
sample=$( awk -v id=$TASK 'NR==id{print $1}' ${INFOFILE} )
mark=$( awk -v id=$TASK 'NR==id{print $2}' ${INFOFILE} )
track=$( awk -v id=$TASK 'NR==id{print $3}' ${INFOFILE} )

export TMP_DIR=$TMP/${sample}${mark}_${RANDOM}
mkdir -p ${TMP_DIR}
echo "[STATUS] Running eval for $sample $mark $track"

cd $ICDIR

# Get and convert the tracks for the blind data:
# ----------------------------------------------
while read -r sample mark track; do
    echo $track
    BLBW=${BLDIR}/blind_${track}
    wget ${baseurl}/blind/${track} -O $BLBW
done < $INFOFILE

# Get chromsizes:
export CHROMFILE=$ANNDIR/hg38.chrom.sizes_main
if [[ ! -s $CHROMFILE ]]; then
    wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/chromInfo.txt.gz -O $CHROMFILE\.gz
    gunzip -c $CHROMFILE.gz | awk -F"\t" -vOFS="\t" '$1 ~ /chr[0-9XY]+$/{print $1,$2}' > $CHROMFILE
fi

# Get and convert the tracks for the blind data:
# ----------------------------------------------
export BLBW=${BLDIR}/blind_${track}
export BLWIG=${BLDIR}/blind_${track%.*}.wig
if [[ ! -s $BLBW ]]; then
    wget ${baseurl}/blind/${track} -O $BLBW

fi

if [[ ! -s $WIGDIR/chr21_blind_${track%%.bigwig}.bedgraph.wig.gz ]]; then
    # Convert to bdg:
    bigWigToBedGraph $BLBW ${BLWIG%%.wig}.bedgraph

    # Convert to chr:
    echo "$sample\t$mark\tblind_${track%%.bigwig}.bedgraph" > $TMP_DIR/convert_info.tsv
    cmd="$CONVBASE ${BLDIR}/ $TMP_DIR/convert_info.tsv $CHROMFILE ${WIGDIR}/"
    echo "$cmd"
    bash -c "${cmd}"
fi


# Get and convert all tracks for the predictions:
# -----------------------------------------------
while read -r num team; do 
    echo "$team -- $num"
    mkdir -p ${R2DIR}/${num}/
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

    export R2BW=${R2DIR}/${num}/${num}_${track}
    if [[ ! -s $R2BW ]]; then
        if [[ $num -eq 200 ]]; then
            echo "[STATUS] Getting + converting: ${parentdir}/${track}"
            wget ${parentdir}/${track} -O $R2BW --no-check-certificate
        else
            # Get the bigwig file:
            echo "[STATUS] Getting + converting: ${parentdir}/${track}"
            wget ${parentdir}/${track} -O $R2BW
        fi
    else
        echo "[STATUS] Bigwig ${track} already exists for team ${team} (${num})"
        ls -sh $R2BW
    fi

    if [[ ! -s $WIGDIR/chr21_${num}_${track%%.*}.bedgraph.wig.gz ]]; then
        export R2BDG=${R2DIR}/${num}/${num}_${track%%.*}.bedgraph
        # Convert to bdg:
        bigWigToBedGraph $R2BW $R2BDG 

        # Convert to chr for Eval: (TODO: Check completion)
        echo "$sample\t$mark\t${num}_${track%%.*}.bedgraph" > $TMP_DIR/convert_info.tsv
        cmd="$CIBASE Convert ${R2DIR}/${num}/ $TMP_DIR/convert_info.tsv $CHROMFILE ${WIGDIR}/"
        echo "$cmd"
        bash -c "${cmd}"
    else
        echo "[STATUS] Converted ${track%%.*} wig already exists for team ${team} (${num})"
        ls -sh $WIGDIR/chr21_${num}_${track%%.*}.bedgraph.wig.gz
    fi

    # Compare the imputed files with the blind data:
    # Vanilla eval:
    # ----------------------------------------------
    export BASEWIG=${num}_${track%%.*}.bedgraph.wig.gz
    echo "$team -- $num -- $BASEWIG"
    if [[ ! -s $EVALFILE ]]; then
        EVALFILE=${STATDIR}/${sample}${mark}_${num}_eval.tsv
        cmd="$CIBASE Eval $WIGDIR blind_${sample}${mark}.bedgraph.wig.gz ${WIGDIR} $BASEWIG $CHROMFILE > ${EVALFILE}"
        echo "$cmd"
        bash -c "$cmd"
    fi
    cat $EVALFILE

    # TODO: Eval with blacklist:

done < $TEAMFILE


# Clean up all tracks:

