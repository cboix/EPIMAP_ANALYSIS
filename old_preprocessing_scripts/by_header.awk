#!/bin/awk
BEGIN { split(cols,out,",")}
NR==1 { for (i=1; i<=NF; i++) ix[$i] = i}
NR>1 { for (i=1; i<length(out); i++){ printf("%s%s", $ix[out[i]], OFS)} print "" }
