#!/usr/bin/R
# -----------------------------------
# Choose a set of random 25bp regions
# -----------------------------------
library(plyr)
domain = system("hostname -d", intern=TRUE)
options(scipen=45)

NREGIONS=2000
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need evalfile and qcfile filenames.")
} else {        
    chromfile = args[1]
    outfile = args[2]
    NREGIONS = as.numeric(args[3])
}

# Read in chromsizes
cs = read.delim(chromfile, sep="\t", header=F)
names(cs) = c('chr','size')

# Concatenate + choose nregions:
cs$nbin = cs$size / 25
totnbin = sum(cs$nbin)

# Sample:
set.seed(1)
reg = as.integer(sample.int(totnbin, size=NREGIONS))
finalbin = c(0, as.integer(cumsum(cs$nbin)))
chrlist = as.character(cs$chr)

# Convert regions back to chromosomes 
outdf = ldply(reg, function(x){
                        s = x - finalbin
                        i = tail(which(s > 0),1)
                        bin = s[i]
                        return(data.frame(chr = chrlist[i], bin=bin)) })

# Write out table:
write.table(outdf, file=outfile, sep="\t", col.names=F, row.names=F, quote=F)
