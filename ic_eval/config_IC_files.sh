#!/bin/bash
# ------------------------------------------------------------
# Setup directories + data for evaluating imputation challenge
# using the relative metrics set
# ------------------------------------------------------------
# Global config and directories:
source ${HOME}/data/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

export ICDIR=${CEDIR}/imputation_challenge
export BLDIR=${ICDIR}/blind
export R2DIR=${ICDIR}/round2
export STATDIR=${ICDIR}/stats
export WIGDIR=${ICDIR}/wig
export OUTDIR=${DBDIR}/out/IC
mkdir -p ${ICDIR} ${BLDIR} ${R2DIR} ${STATDIR} ${WIGDIR} ${OUTDIR}

# Files for run:
export baseurl="http://mitra.stanford.edu/kundaje/ic"
export ciurl="https://public.hoffman2.idre.ucla.edu/ernst/C86ZP"
export qnormurl="http://mitra.stanford.edu/kundaje/jmschr/proj/2020_encode_imputation_challenge/qnorm_bigwigs/"
export CIBASE="java -mx8000M -jar $BINDIR/ChromImpute_SRC/ChromImpute_coeffv.jar"
export INFOFILE=${ICDIR}/ic_blind_tracks.tsv
export TEAMFILE=${ICDIR}/team_name_round2.tsv

# Get team information:
if [[ ! -s $TEAMFILE ]]; then
    wget ${baseurl}/round2/team_name_round2.tsv -O $TEAMFILE
    echo "200\tChromImpute" >> $TEAMFILE
fi

# Get chromsizes:
export CHROMFILE=$ANNDIR/hg38.chrom.sizes_main
if [[ ! -s $CHROMFILE ]]; then
    wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/chromInfo.txt.gz -O $CHROMFILE\.gz
    gunzip -c $CHROMFILE.gz | awk -F"\t" -vOFS="\t" '$1 ~ /chr[0-9XY]+$/{print $1,$2}' > $CHROMFILE
fi

# source $BINDIR/ic_eval/process_IC_files.sh

