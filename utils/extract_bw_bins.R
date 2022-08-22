#!/usr/bin/R
# ----------------------------------
# Extract a set of bins from a file:
# ----------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    source('~/data/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
} else {
    source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
}
library(dplyr)
options(scipen=45)

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need evalfile and qcfile filenames.")
} else {        
    selfile = args[1]
    infofile = args[2]
    regionsfile = args[3]
    outprefix = args[4]
}

chunk = as.integer(as.numeric(Sys.getenv("SGE_TASK_ID")))

# Read in sample tables:
seldf = read.delim(selfile, sep="\t", header=F, stringsAsFactors=F)
names(seldf) = c('id','mark')

# Choose set:
mark = as.character(seldf$mark[chunk])
id = as.character(seldf$id[chunk])
print(paste(mark, id))

# Get prefixes:
sampdf = read.delim(infofile, sep="\t", header=F)
names(sampdf) = c('id', 'mark', 'obsfile')
idx = which(sampdf$id == id & sampdf$mark == mark)
obs.pref = sampdf$obsfile[idx[1]]
print(obs.pref)

#obs.pref = filtsamp$obsfile[1]
imp.pref = paste0('impute_', id, '_', mark)
obsdir = 'ChromImpute/converted/'
impdir = 'ChromImpute/imputed/'

# Read in random regions 
rdf = read.delim(regionsfile, sep="\t", header=F, stringsAsFactors=F)
names(rdf) = c('chr','bin')
chrlist = sort(unique(as.character(rdf$chr)))

# Functions:
read.track = function(dir, chr, pref, ind){
    filename = paste0(dir, chr, '_', pref, '.wig.gz')
    as.numeric(scan(gzfile(filename), 'c', skip=2))[ind]
}

# Output files
impfile = paste0(outprefix, '_imp_', chunk, '.tsv')
obsfile = paste0(outprefix, '_obs_', chunk, '.tsv')

# Decide whether to (re)run:
get.regions = FALSE
if (file.exists(impfile)){
    x = scan(impfile, 'c', quiet=T)
    if (length(x) != nrow(rdf)){ get.regions = TRUE }
} else { get.regions = TRUE }

if (file.exists(obsfile)){
    x = scan(obsfile, 'c', quiet=T)
    if (length(x) != nrow(rdf)){ get.regions = TRUE }
} else { get.regions = TRUE }


# Read all regions + write them:
if (get.regions){
    for (chrom in chrlist){
        bool.append = (chrom != chrlist[1])
        print(paste('Reading in indices in chromosome:', chrom))
        subdf = filter(rdf, chr == chrom)
        ind = as.integer(as.character(subdf$bin))

        # Read imputed bins
        imp = read.track(impdir, chrom, imp.pref, ind)
        # Write imputed bins
        write.table(imp, file=impfile, quote=F, sep="\t", row.names=F,
                    col.names=F, append=bool.append)

        # Read observed bins
        obs = read.track(obsdir, chrom, obs.pref, ind)
        # Write observed bins
        write.table(obs, file=obsfile, quote=F, sep="\t", row.names=F,
                    col.names=F, append=bool.append)
    }
}
