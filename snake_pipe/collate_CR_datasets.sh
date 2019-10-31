#!/bin/bash
# Combine data to make datasets for chromatin regulators
export PROJ=$1 # CRs or EPIMAP_ANALYSIS
export DOMAIN=$(hostname -d)
if [[ $DOMAIN == 'broadinstitute.org' ]]
then
    export PROJDIR=$HOME/data/${PROJ}/db
    export DBDIR=$HOME/data/EPIMAP_ANALYSIS/db
    export BINDIR=$HOME/data/EPIMAP_ANALYSIS/bin
else
    export PROJDIR=$HOME/${PROJ}/db
    export DBDIR=$HOME/EPIMAP_ANALYSIS/db
    export BINDIR=$HOME/EPIMAP_ANALYSIS/bin
fi
export ANNDIR=$DBDIR/Annotation
export CHDIR=$PROJDIR/ChromHMM/binarized
export BDGDIR=$PROJDIR/ChromImpute/converted
export TMPDIR=${TMP}/${RANDOM}  # Work in TMPDIR
mkdir -p $CHDIR/full $BDGDIR/full $TMPDIR $TMPDIR/full

# ===========
# CR datasets
# ===========
cd $CHDIR
if [[ "$PROJ" == "CRs" ]]
then
    ls -d */*_chr1_binary.txt.gz | awk '{sub("_chr1_binary.txt.gz","",$0); print $0}' | sort -u > $CHDIR/prefix_list
    awk '{sub(".*/","",$0); print $0}' $CHDIR/prefix_list | sort -u > $CHDIR/cell_types
    awk '{sub(".*/BSS","BSS"); sub("_.*",""); print $0}' $CHDIR/prefix_list | sort -u > $CHDIR/reduced_cell_types
    # echo "A549\nGM12878\nH1-hESC\nHepG2\nK562" > $CHDIR/reduced_cell_types
elif [[ "$PROJ" == "EPIMAP_ANALYSIS" ]]
then
    ls -d */*_full_binary.txt.gz | awk '{sub("_full_binary.txt.gz","",$0); print $0}' | sort -u > $CHDIR/prefix_list
    awk '{sub(".*/","",$0); print $0}' $CHDIR/prefix_list | sort -u > $CHDIR/cell_types
    awk '{sub(".*/BSS","BSS"); sub("_.*",""); print $0}' $CHDIR/prefix_list | sort -u > $CHDIR/reduced_cell_types
    # echo "A549\nGM12878\nH1-hESC\nHepG2\nK562" > $CHDIR/reduced_cell_types
fi

while read celltype; do
    echo "Binned data for ${celltype}"
    ENDINGS=( "${celltype}_c_t" "${celltype}_c.sub_t.sub" "${celltype}_c_t.sub" )
    for CTEXT in ${ENDINGS[@]}; do
        while read chr; do
            FINAL_CHRFILE=${CHDIR}/full/${CTEXT}_${chr}_CRs_binary.txt
            if [[ ! -s ${FINAL_CHRFILE}.gz ]]; then
                echo "-Making files for chromosome $chr for $CTEXT"
                while read prefix; do
                    echo $prefix
                    CHRFILE=${TMPDIR}/full/${CTEXT}_${chr}_CRs_binary.txt
                    zcat "${prefix}_${chr}_binary.txt.gz" | pr -mts ${CHRFILE} - > ${CHRFILE}\.tmp
                    mv ${CHRFILE}\.tmp ${CHRFILE}
                done < <( grep "/${CTEXT}$" $CHDIR/prefix_list)
                # Replace old with new copy:
                mv ${CHRFILE} ${FINAL_CHRFILE}
                gzip ${FINAL_CHRFILE} -f
            fi 
        done < <( awk '{print $1}' ${ANNDIR}/hg19.chrom.sizes_main )
    done
done < $CHDIR/reduced_cell_types

# ONLY HISTONE MARKS:
while read celltype; do
    echo "Binned data for ${celltype}"
    ENDINGS=( "${celltype}_c_t" "${celltype}_c.sub_t.sub" "${celltype}_c_t.sub" )
    for CTEXT in ${ENDINGS[@]}; do
        while read chr; do
            FINAL_CHRFILE=${CHDIR}/full/${CTEXT}_${chr}_marks_binary.txt
            if [[ ! -s ${FINAL_CHRFILE}.gz ]]; then
                echo "-Making files for chromosome $chr for $CTEXT"
                while read prefix; do
                    echo $prefix
                    CHRFILE=${TMPDIR}/full/${CTEXT}_${chr}_marks_binary.txt
                    zcat "${prefix}_${chr}_binary.txt.gz" | pr -mts ${CHRFILE} - > ${CHRFILE}\.tmp
                    mv ${CHRFILE}\.tmp ${CHRFILE}
                done < <( grep "/${CTEXT}$" $CHDIR/prefix_list | awk '$1 ~ /^H[0-4]/ || $1 ~ /DNase/' )
                # Replace old with new copy:
                mv ${CHRFILE} ${FINAL_CHRFILE}
                gzip ${FINAL_CHRFILE} -f
            fi 
        done < <( awk '{print $1}' ${ANNDIR}/hg19.chrom.sizes_main )
    done
done < $CHDIR/reduced_cell_types


# -------------------------
# Bedgraphs (25/50/100/200)
# CellTypes (all combinations)
# -------------------------
# Both FC and log.10.Pval:
# NOTE if CRs, go to BDGDIR
cd $BDGDIR
if [[ "$PROJ" == "CRs" ]]
then 
    ls -d */chr1_*.pval.signal.bedgraph.gz.wig.gz | awk '{sub("/chr1_FINAL_","\t",$0); sub("_VS_","\t",$0); sub($1"_","",$2); sub("FINAL_WCE_","\t",$0); sub(".pval.*$","\t",$0); print $0}' | sort -u > $BDGDIR/prefix_list
else
    ls -d chr1_*.pval.signal.bedgraph.gz.wig.gz | awk '{sub("chr1_FINAL_","",$0); sub(".pval.*$","",$0); sub("_VS_FINAL_WCE_","\t",$0); sub("_"$2,"\t"$2,$0); print $0,ct}' | sort -u > $BDGDIR/prefix_list
fi

CT=( "A549" "GM12878" "H1-hESC" "HepG2" "K562" )
SIGEND="signal.bedgraph.gz.wig.gz"
# rm $BDGDIR/full/*.wig.gz -rf
for CELLTYPE in ${CT[@]}
do
    echo "Collating bedgraphs for ${CELLTYPE}"
    # chr="chr21"
    while read chr
    do
        while read ept TAcell CTcell
        do
            echo "Bedgraph data for ${ept} ${TAcell} vs. ${CTcell}"
            if [[ "$PROJ" == "CRs" ]]
            then 
                FC=${BDGDIR}/${ept}/${chr}_FINAL_${ept}_${TAcell}_VS_FINAL_WCE_${CTcell}.fc.${SIGEND}
                PVAL=${BDGDIR}/${ept}/${chr}_FINAL_${ept}_${TAcell}_VS_FINAL_WCE_${CTcell}.pval.${SIGEND}
            else
                FC=${BDGDIR}/${chr}_FINAL_${ept}_${TAcell}_VS_FINAL_WCE_${CTcell}.fc.${SIGEND}
                PVAL=${BDGDIR}/${chr}_FINAL_${ept}_${TAcell}_VS_FINAL_WCE_${CTcell}.pval.${SIGEND}
            fi
            FCFULL=${TMPDIR}/${TAcell}_VS_${CTcell}_${chr}_CRs.fc.25.wig
            PVALFULL=${TMPDIR}/${TAcell}_VS_${CTcell}_${chr}_CRs.pval.25.wig
            zcat ${FC} | pr -mts ${FCFULL} - > ${FCFULL}\.tmp
            mv ${FCFULL}\.tmp ${FCFULL}
            zcat ${PVAL} | pr -mts ${PVALFULL} - > ${PVALFULL}\.tmp
            mv ${PVALFULL}\.tmp ${PVALFULL}
        done < <( grep "$CELLTYPE" $BDGDIR/prefix_list )
    done < <( awk '{print $1}' ${ANNDIR}/hg19.chrom.sizes_main )
done

# Move files to real directory:
cd $TMPDIR
while read file
do
    echo "Compressing $file"
    # Replace old with new (compressed):
    awk 'NR==1{gsub("track type=wiggle_0 name=","",$0); gsub("_observed","",$0); print $0}NR==2{sub("\tfixedStep.*$","",$0); print $0}NR > 2{print $0}' ${TMPDIR}/${file} | gzip -c > ${BDGDIR}/full/${file}.gz
    rm ${TMPDIR}/${file}
done < <( cd $TMPDIR; ls *_CRs.*.wig )


rm -rf ${TMPDIR}
