#!/bin/bash
# For running snakemake pipeline:
# Does: 
# Looks at all outfiles
# - Name
# - Progress percent if finished.
# - Checks common errors (memory, MissingInput, bam sorting).
# - moves all annotated files to directories. (after putting them in a list)
# So that we can check the others for issues...

cd $DBDIR/out/ChIP
while read file
do
    if [[ "1" == "$( grep "Nothing to be done" $file -c )" ]]
    then
        echo "$file was done"
        rm $file
    else
        if [[ "1" == "$( grep "(100%) done" $file -c )" ]]
        then
            echo "$file now finished"
            rm $file
        else
            if [[ "1" == "$( grep "OSError: Cannot call rmtree" $file -c )" ]]
            then
                echo "$file bad symlink?"
                echo $file >> cleanup_files
            fi
        fi
    fi
done < <( ls *.o* )


# For ChromImpute pipeline:
while read file
do
    if [[ "1" == "$( grep "Finished DISTANCE sucessfully" $file -c )"  ]]
    then
        echo "$file is done"
        rm $file
        # if [[ "4" == "$( wc -l $file | awk '{print $1}' )" ]]
        # then
        #     echo $file > tasks_run
        # fi
    fi
done < <( ls dist*.o* | sort -n )

# For ChromImpute pipeline:
while read file
do
    if [[ "1" == "$( grep "Finished FEATURES sucessfully" $file -c )"  ]]
    then
        echo "$file is done"
        rm $file
        # if [[ "4" == "$( wc -l $file | awk '{print $1}' )" ]]
        # then
        #     echo $file > tasks_run
        # fi
    fi
done < <( ls feat*.o* | sort -n )

# RERUN ALL FAILED DISTANCES....
while read TASK
do
    qsub -cwd -l h_vmem=10G -N dist_${TASK}_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh DISTANCE $TASK"
done < tasks_run


# For ChromImpute pipeline:
cd ${DISTANCE_DIR}
while read file
do
    if [[ "0" == "$( wc -l $file | awk '{print $1}' )"  ]]
    then
        f2=${file##*_}
        MARK=${f2%.txt}
        SAMPLE=${file%_*}
        TASK=$( grep -n "${MARK}_${SAMPLE}.sub" $SAMPLEMARK_TAB | awk -vFS=":" '{print $1}' )
        echo "$TASK : $SAMPLE with $MARK"

        # SUBMIT IT:
        qsub -cwd -l h_vmem=10G -N dist_${TASK}_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh DISTANCE $TASK"
    fi
done < <( ls *.txt | sort -n )


while read link 
do
    echo $link
    wget ${link} 
done < TX_links.tsv

# TODO MAKE THIS BETTER.
# RERUN ALL FAILED PREDICTORS.
while read TASK
do
    qsub -cwd -l h_vmem=15G -q long -N pred_${TASK}_${UQID} -j y -b y -V -r y -o $DBDIR/out/ChromImpute "$BINDIR/run_ChromImpute.sh PREDICTORS $TASK"
done < task_run

