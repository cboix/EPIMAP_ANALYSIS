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

grep "Nothing to be done" -l *.o* > previously_finished
grep "(100%) done" -l *.o* > finished_runs
# grep "%) done" testfile | awk '{sub(".*\\(",""); sub("%\\) done.*",""); print $0}'
grep "(100%) done" -l *.o* > finished_runs
grep "OSError: Cannot call rmtree" -l *.o* > tree_oserror
grep "OSError:" -l *.o* > other_oserror
grep "No space " -l *.o* > memory_error
grep "FileNotFound" -l *.o* > missingInput_files
grep "unsorted positions" -l *.o* > sort_error
grep "unsorted positions" -l *.o* | sort +0 -1 | awk 'BEGIN{print "name","sort_error"}{print $0,1}' > sort_error

# *** Percent status of files *** 
while read file 
do 
    tac $file | grep -m1 "%) done" - | awk -v file=$file -vOFS="\t" '{sub(".*\\(",""); sub("%\\) done.*",""); print file,$0}' >> pct_status.log
done < <( ls *.o*) 

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
