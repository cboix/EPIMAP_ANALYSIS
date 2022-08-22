#!/bin/bash
# ==========================================================
# Get IHEC bigWig files and turn into bedgraphs for ChromHMM
# Gets BW to IHEC dir, and then redirect bedgraphs to the correct ChIP directories
# ==========================================================
donid=$( sed "${SGE_TASK_ID}q;d" ${ANNDIR}/IHEC_donors.tsv )
echo "Processing Donor: $donid"
INFO=$ANNDIR/IHEC_bw_histonemods.tsv

IFS=$'\t'
while read link consortium assembly file sampleid donor factor processing
do
    # Put files into consortium specific directories
    CONSDIR=${IHECDIR}/${consortium}
    mkdir -p $CONSDIR

    # Call input "WCE"
    if [[ "${factor}" == "Input" ]]
    then
        factor="WCE"
    fi
    BWFILE=${CONSDIR}/${sampleid}_${donor}_${factor}_${assembly}_${processing}.bigWig

    BEDFILE=${CONSDIR}/${sampleid}_${donor}_${factor}_${assembly}_${processing}.bedgraph

    if [[ ! -s ${BEDFILE}.gz ]]
    then
        if [[ ! -s ${BWFILE} ]]
        then
            echo "Downloading ${BWFILE}"
            wget ${link} -o ${BWFILE}.log -O ${BWFILE}
        fi
        echo "Making ${BEDFILE}"
        bigWigToBedGraph ${BWFILE} ${BEDFILE}

        # Trim excess precision that doesnt really exist!
        awk '{printf("%s\t%s\t%s\t%.0f\n",$1,$2,$3,$4)}' ${BEDFILE} | awk -vOFS="\t" 'NR==1{chr=$1; start = $2; end=$3; a = $4}NR > 1{ if( end == $2 && a == $4 ){ end = $3; } else { print chr,start,end,a; chr = $1; start = $2; end = $3; a = $4; }} END{print chr,start,end,a}' | gzip -c > ${BEDFILE}.gz
        rm -f ${BEDFILE} ${BWFILE}
    fi
done < <( awk -vFS="\t" -v don=$donid '$6 == don' $INFO )


# TODO Put bedfile in correct directory
# TODO Map hg19 to hg38
