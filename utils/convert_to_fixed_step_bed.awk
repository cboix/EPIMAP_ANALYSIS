#!/bin/awk
# -----------------------------------------------------------
# Converts a bedgraph into a fixed-step bed
# NOTE: Need to further process to make a fixed-step wig file
# -----------------------------------------------------------
# TODO: Currently doesn't avoid going over chrom-size 
# Needs to read in secondary file into list before.

BEGIN {
    chr="chr0"; 
    s=1; v=0; l=0;
}

$1 == chrom {
    chr1=$1; s1=$2+1; e1=$3; v1=$4; l1=e1-s1+1;
    # If any remnant from the last line, amalgamate:
    if (chr == chr1){
        if (l > 0){
            lmiss = step - l
            if (l1 >= lmiss){
                # Print the average if enough for full step:
                avg = (v * l + v1 * lmiss) / step;
                print chr1, s, s + step - 1, avg;
                l1 = l1 - lmiss; v1 = v1; s1 = s + step;
            } else {
                # Carry value over if still falls short:
                v1 = (v * l + v1 * l1) / (l + l1);
                l1 = l + l1;
                s1 = s;
            }
        }
    }
    s = s1; l = l1; v = v1; chr = chr1;
    # Shear off any full steps if possible:
    while(l >= step){
        s2 = s + step;
        print chr, s, s2 - 1, v;
        l = e1 - s2+1;
        s = s2;
        s1 = s2;
    }
}


