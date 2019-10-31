#!/bin/bash
# ----------------------------------------------
# Inner command for enrichments (modules or raw)
# ----------------------------------------------
echo "[STATUS] Running inner cmd: $COMMAND"

# Make temporary directories:
export TMP_DIR=${TMP}/CTOM_${COMMAND}_${MODEL}_${RANDOM}
mkdir -p ${TMP_DIR}

# Make infofile for great/fimo + others?:
if [[ "$COMMAND" == "GREAT" ]] || [[ "$COMMAND" == "FIMO" ]] || [[ "$COMMAND" == "GWAS" ]]; then
    echo "[STATUS] Making infofile for command $COMMAND"
    INFOFILE=${TCDIR}/enrichment_infofile.tsv
    ls $SPLITPREFIX*.bed | awk -vOFS="\t" -v sp=$SPLITPREFIX '{a=$1; sub(".bed","",a); sub(sp,"",a); print $1, a}' > ${INFOFILE}
    head ${INFOFILE}
    wc -l ${INFOFILE}
    INUM=$( cat $INFOFILE | wc -l )
fi

# Clustering and extraction variables/arguments:
if [[ "$COMMAND" == "EXTRACT" ]] || [[ "$COMMAND" == "CLUSTER" ]]; then
    echo "[STATUS] Setting arguments for command $COMMAND"
    CPPREF=${EOUTPREF}_r200_e0_allchr
    # IF SMERGE=2, use mergedfile. (+ _merged)
    if [[ "$SMERGE" == "2" ]]; then 
        CPPREF=${CPPREF}_merged
    fi
    CPFILE=${CPPREF}_csr.cp.gz
    ATTRFILE=${CPPREF}_attr.cp.gz
    MAINARGS="--filename ${CPFILE} --attrfile $ATTRFILE --outprefix ${COUTPREF} --K ${NCLUST}"
    # Add option of merging or not:
    if [[ "$SMERGE" != "0" ]]; then
        MAINARGS="$MAINARGS --mergestates"
        if [[ "$SMERGE" == "1" ]]; then
            MAINARGS="$MAINARGS --mergeset $MERGED_STATES"
        fi
    fi
    # Add secondary (raw signal) file:
    if [[ "$MWITH" == "1" ]]; then
        H3DIR=${IMPUTED_DIR}/H3K27ac
        HPREFIX=H3K27ac_all_bin_${JOBSUF}
        H3PREF=${H3DIR}/$HPREFIX
        H3_CPPREF=${H3PREF}_r25_e${EXTSIZE}_allchr_merged
        H3_CPFILE=${H3_CPPREF}_csr.cp.gz
        H3_ATTRFILE=${H3_CPPREF}_attr.cp.gz
        MAINARGS="$MAINARGS --h3file ${H3_CPFILE} --h3attr ${H3_ATTRFILE}"
    fi
    # Add keep idx if previously specified (in ipynb):
    if [[ "$ELEMENT" == "ENH" ]]; then
        MAINARGS="$MAINARGS --keepfile $DML_DIR/ENH_masterlist_indices_ind.cp.gz"
    elif [[ "$ELEMENT" == "PROM" ]]; then
        MAINARGS="$MAINARGS --keepfile $DML_DIR/PROM_masterlist_indices_ind.cp.gz"
    fi
fi

if [[ "$COMMAND" == "PREPROCESS" ]]; then
    SAMPLE=$( cut -f1 ${SAMPLEMARK_TAB} | sort -u | sed "${TASK}q;d" - | awk -v FS="\t" '{print $1}' )
    MAINCMD="source activate pytorch_env > /dev/null; python ${BINDIR}/aggregate_binary_txt_to_npy.py preprocess --prefix $SAMPLE $MAINARGS"
    echo "[STATUS] Running $COMMAND command:"
    echo "$MAINCMD"
    bash -c "$MAINCMD"

elif [[ "$COMMAND" == "AGGREGATE" ]]; then
    MAINCMD="source activate pytorch_env > /dev/null; python ${BINDIR}/aggregate_binary_txt_to_npy.py main $MAINARGS"
    echo "[STATUS] Running $COMMAND command:"
    echo "$MAINCMD"
    bash -c "$MAINCMD"

elif [[ "$COMMAND" == "MAKEMAT" ]]; then
    echo "[STATUS] Making matrix for states $STATES (${ELEMENT}) for model ${MODEL}"
    # Secondary (raw signal) file:
    H3DIR=${IMPUTED_DIR}/H3K27ac
    HPREFIX=H3K27ac_all_bin_${JOBSUF}
    H3PREF=${H3DIR}/$HPREFIX
    H3_CPPREF=${HPREF}_r25_e${EXTSIZE}_allchr_merged
    H3_CPFILE=${H3_CPPREF}_csr.cp.gz
    H3_ATTRFILE=${H3_CPPREF}_attr.cp.gz

elif [[ "$COMMAND" == "EXTRACT" ]]; then
    echo "[STATUS] EXTRACTING elements for ${ELEMENT} for model ${MODEL}"
    # Main command + extract:
    MAINCMD="source activate pytorch_env > /dev/null; python ${BINDIR}/cluster_binary_mat.py main $MAINARGS --extract"
    echo "[STATUS] Running $COMMAND command:"
    echo "$MAINCMD"
    bash -c "$MAINCMD"

    # Split files for later running enrichments:
    # Modified command relative to clustering:
    NBED=$(ls $TCDIR/raw_*_assignments.bed | wc -l )
    NCOUNTED=$( wc ${COUTPREF}_counts.tsv | awk '{print $1}' )
    echo "[STATUS] $NBED BEDfiles, and $NCOUNTED counted"
    if [[ "$NBED" != "$NCOUNTED" ]]; then
        echo "[STATUS] Separating files for ${FPREFIX} - ${CPREFIX}"
        SEPPREFIX=$MTDIR/${FPREFIX}
        # Ensure correct bedfiles matching DHS locations
        COORDMAP=${DML_DIR}/masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r200_e0_names.core.srt.tsv
        COORDSRTMAP=${DML_DIR}/masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r200_e0_names.core.srt.tsv.reord
        awk -vOFS="\t" '{print $4,$1,$2,$3}' $CORECOORD | sort +0 -1 > $TMP_DIR/corecoord_srt.tsv
        if [[ ! -s $COORDSRTMAP ]] || [[ $COORDMAP -nt $COORDSRTMAP ]]; then
            awk -vOFS="\t" '{print $2, $3}' $COORDMAP | sort +0 -1 > $COORDSRTMAP
        fi
        rm ${COUTPREF}_nbp.tsv ${COUTPREF}_counts.tsv
        while read file; do 
            IDNUM=$( awk -F"\t" -vOFS="\t" 'NR == 2{print $1}' $file )
            echo "[STATUS] Converting file $(basename $file) - $IDNUM"
            # Output file:
            SEPFILE=${SEPPREFIX}${IDNUM}.bed
            awk -F"\t" -vOFS="\t" 'NR > 1{print $2, $1}' $file |  sort +0 -1 | join - ${COORDSRTMAP} | awk -vOFS="\t" '{print $3,$2}' | sort +0 -1 | join ${TMP_DIR}/corecoord_srt.tsv - | awk -vOFS="\t" '{print $2,$3,$4,$5}' | sort -k1V -k2n > ${SEPFILE}
            awk -vOFS="\t" -v id=$IDNUM '{a=a+ $3-$2}END{print id, a}' ${SEPFILE} >> ${COUTPREF}_nbp.tsv
            awk -vOFS="\t" -v id=$IDNUM 'END{print id, NR}' ${SEPFILE} >> ${COUTPREF}_counts.tsv
        done < <(ls $TCDIR/raw_*_assignments.bed)
    fi

elif [[ "$COMMAND" == "CLUSTER" ]]; then
    echo "[STATUS] Clustering elements for ${ELEMENT} for model ${MODEL}"
    # Main command for clustering:
    MAINCMD="source activate pytorch_env > /dev/null; python ${BINDIR}/cluster_binary_mat.py main $MAINARGS"
    echo "[STATUS] Running $COMMAND command:"
    echo "$MAINCMD"
    bash -c "$MAINCMD"

    # Split files for later running enrichments:
    if [[ ! -s ${CLSBED%bed}fixed.bed ]] || [[ $CLSBED -nt ${CLSBED%bed}fixed.bed ]]; then
        # Ensure correct bedfiles matching DHS locations
        COORDMAP=${DML_DIR}/masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r200_e0_names.core.srt.tsv
        COORDSRTMAP=${DML_DIR}/masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r200_e0_names.core.srt.tsv.reord
        if [[ ! -s $COORDSRTMAP ]] || [[ $COORDMAP -nt $COORDSRTMAP ]]; then
            awk -vOFS="\t" '{print $2, $3}' $COORDMAP | sort +0 -1 > $COORDSRTMAP
        fi
        awk -F"\t" -vOFS="\t" 'NR > 1{print $2, $1}' $CLSBED |  sort +0 -1 | join - ${COORDSRTMAP} | awk -vOFS="\t" '{print $3,$2}' | sort +0 -1 > $TMP_DIR/tmpcls.bed
        awk -vOFS="\t" '{print $4,$1,$2,$3}' $CORECOORD | sort +0 -1 > $TMP_DIR/corecoord_srt.tsv
        join ${TMP_DIR}/corecoord_srt.tsv ${TMP_DIR}/tmpcls.bed | awk -vOFS="\t" '{print $2,$3,$4,$5}' | sort -k1V -k2n > ${CLSBED%bed}fixed.bed
        # Separate files into directory:
        echo "[STATUS] Separating files for ${FPREFIX} - ${CPREFIX}"
        awk '$1 == "chr17"' ${CLSBED%bed}fixed.bed | head -n 2
        awk -vOFS="\t" -v prefix=${MTDIR} -v fpref=${FPREFIX} '{print $1,$2,$3 > prefix"/"fpref$4".bed"}' ${CLSBED%bed}fixed.bed
        awk -vOFS="\t" 'BEGIN{split("",a,"");}{a[$4]=a[$4] + $3-$2}END{b=length(a) -1; for(i=0; i<=b;i++){print i, a[i]}}' ${COUTPREF}_assignments.fixed.bed > ${COUTPREF}_nbp.tsv
        awk 'NR > 1{print $1}' ${CLSBED} | sort -n | uniq -c | awk -vOFS="\t" '{print $2, $1}' > ${COUTPREF}_counts.tsv
    fi

    # Plot the clusters alone:
    source activate mv_env > /dev/null
    plotcmd="R --slave -f ${BINDIR}/plot_modules_only.R --args ${COUTPREF} '$TAGLINE' $TCIMG"
    echo "[STATUS] Plotting computed clusters"
    echo "$plotcmd"
    bash -c "$plotcmd"
    source deactivate > /dev/null

elif [[ "$COMMAND" == "PLOTCLS" ]]; then
    # Plot the clusters alone:
    source activate mv_env > /dev/null
    plotcmd="R --slave -f ${BINDIR}/plot_modules_only.R --args ${COUTPREF} '$TAGLINE' $TCIMG"
    echo "[STATUS] Plotting computed clusters"
    echo "$plotcmd"
    bash -c "$plotcmd"
    source deactivate > /dev/null


elif [[ "$COMMAND" == "GREAT" ]]; then
    # Submit great pipeline, no background:
    FINALJOB=$( ${BINDIR}/great_enrichment/submit_GREAT_pipeline.sh ${INFOFILE} ${GTDIR} NONE | awk '{a=$1}END{print $1}')
    echo "[STATUS] Final job in pipe was $FINALJOB"

    # Plot above results + modules
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=3:00:00 -N plot_great_${UQID} -hold_jid ${FINALJOB} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env >> /dev/null; R --slave -f ${BINDIR}/plot_modules_with_great.R --args ${GTDIR} ${COUTPREF} '$TAGLINE' $TCIMG; source deactivate >> /dev/null"

    # Submit great command, with background:
    # BKGFILE=${COUTPREF}_assignments.fixed.bed  # All enhancers in any cluster
    # qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=3:00:00 -N great_${FPREFIX}_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "${BINDIR}/great_enrichment/submit_GREAT_pipeline.sh ${INFOFILE} ${GTBKDIR} ${BKGFILE}"

elif [[ "$COMMAND" == "MOTIF" ]]; then
    # Use above output to run the motif pipeline:
    mkdir -p $MTDIR/results
    ENRFILE=${MTDIR}/results/motif_enrichment/merged_enrichments_0_-99.tsv
    RUNMOT=1
    if [[ ! -s $ENRFILE ]] || [[ "$RUNMOT" == "1" ]]; then
        echo "[STATUS] Queueing motif pipeline"
        FINALJOB=$( $MVDIR/submit_enrichment_motifs.sh -i ${MTDIR}/ -o "hg19" -d $MTDIR/results/ | awk '{a = $1}END{print a}')
    else
        FINALJOB="null_name"
    fi
    echo "[STATUS] Final job id is: $FINALJOB"

    # Plot enrichment results:
    echo "[STATUS] Queueing plot motifs with centers"
    qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=3:00:00 -N plot_motifs_${UQID} -hold_jid ${FINALJOB} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env >> /dev/null; R --slave -f ${BINDIR}/plot_modules_with_motifs.R --args ${ENRFILE} ${COUTPREF} '$TAGLINE' $TCIMG; source deactivate >> /dev/null"


elif [[ "$COMMAND" == "FIMO" ]]; then
    # Submit great pipeline, no background:
    FINALJOB=$( ${BINDIR}/meme_enrichment/submit_fimo_enrichment.sh ${INFOFILE} ${MTDIR} | awk '{a=$1}END{print $1}')
    echo "[STATUS] Final job in pipe was $FINALJOB"

    # Plot above results + modules
    # qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=3:00:00 -N plot_great_${UQID} -hold_jid ${FINALJOB} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env >> /dev/null; R --slave -f ${BINDIR}/plot_modules_with_great.R --args ${GTDIR} ${COUTPREF} '$TAGLINE' $TCIMG; source deactivate >> /dev/null"

elif [[ "$COMMAND" == "GWAS" ]]; then
    EXTLIST=(0 1000 5000)
    for extgwas in ${EXTLIST[@]}; do
        echo $extgwas
        JOBTAG=gwas_${FPREFIX}_${extgwas}_${UQID}
        GWFILE=${GWPREFIX}_${extgwas}_enrich.tsv

        # Run enrichment calculation:
        if [[ ! -s $GWFILE ]]; then
            # NOTE: 3 hrs for 300, queue for 8 or so for 800+ ? 
            echo "[STATUS] Queuing gwas enrichment"
            qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=8:00:00 -N ${JOBTAG} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "java -mx8000M -Djava.io.tmpdir=${TMP_DIR} -jar ${GWASJAR} $SPLITPREFIX $GWPREFIX $GWASCATALOG $extgwas $INUM"
        fi

        # Plot enrichment:
        echo "[STATUS] Queuing plotting gwas enrichment with clusters"
        qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:00:00 -N plot_${JOBTAG} -hold_jid ${JOBTAG} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env >> /dev/null; R --slave -f ${BINDIR}/plot_modules_with_gwas.R --args ${GWFILE} ${COUTPREF} '$TAGLINE' $extgwas $TCIMG; source deactivate >> /dev/null"

        echo "[STATUS] Queuing plotting gwas enrichment alone"
        qsub -cwd -P compbio_lab -l h_vmem=10G -l h_rt=1:00:00 -N plot_${JOBTAG} -hold_jid ${JOBTAG} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env >> /dev/null; R --slave -f ${BINDIR}/gwas_enrichment/plot_gwas_alond.R --args ${GWFILE} ${COUTPREF} '$TAGLINE' $extgwas $TCIMG; source deactivate >> /dev/null"
    done

elif [[ "$COMMAND" == "LDAK" ]]; then
    break

elif [[ "$COMMAND" == "LDSC" ]]; then
    break

else
    echo "[STATUS] Unknown chmm to modules command: $COMMAND"
    echo "Please use one of: PREPROCESS, AGGREGATE, MAKEMAT, CLUSTER."
    echo "or for enrichment: GREAT, MOTIF, GWAS, LDAK, LDSC"
    exit
fi

cd $TCDIR
tar -cvzf $TCDIR/images.tar.gz img
echo "[STATUS] Images at $TCDIR/images.tar.gz"

rm -rf ${TMP_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Finished $COMMAND sucessfully in $runtime seconds."
