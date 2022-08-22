#!/usr/bin/R
# ------------------------------------------------------
# Create the permuted catalogs for significance analysis
# ------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))

# Arguments for loading data:
usetree = 'enhancers'
# usetree = 'roadmap'  # For old epi only
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only) }
source(paste0(bindir, 'load_gwastree_analysis.R'))

# Under dbdir:
gtdir = "gwas_tree_analysis/"
perdir = paste0(gtdir, "permuted_catalogs/")
cmd = paste('mkdir -p', perdir)
system(cmd)

# ----------------------------------------
# Create the permuted catalogs
# ----------------------------------------
gwcatfile = 'gwascatalog_may03_2019_noquotes.txt'
gwrdafile = sub(".txt", ".Rda", gwcatfile)
gwintrdafile = sub(".txt", "_intersections.Rda", gwcatfile)

load(gwrdafile)
load(gwintrdafile)
nhitdf = aggregate(subjectHits ~ uid, qdf, length)
nhitdf = aggregate(uid ~ subjectHits, nhitdf, length)
head(nhitdf)

# Separate by number of SNPs, not intersections:
nsnpdf = aggregate(pValue ~ uid, gwdf, length)
nsnpdf = aggregate(uid ~ pValue, nsnpdf, length)
TOTSNP = nrow(gwdf)
head(nsnpdf)
NCAT = 1000

for (j in 1:nrow(nsnpdf)){
    # Separate by number of SNPs, not intersections:
    NSNPS = nsnpdf$pValue[j] 
    # Concatenate NCAT permuted intersection sets:
    gwperfile = paste0(perdir, sub(".txt", paste0("_", NCAT, sprintf("_%06d.Rda", NSNPS)), gwcatfile))
    if (!file.exists(gwperfile)){
        # TODO: WRONG - FIXME should be by NUMBER OF SNPs:
        print(NSNPS)
        # Select permutation sets:
        set.seed(1414)
        hits = sapply(rep(NSNPS, NCAT), function(x){sample(1:TOTSNP, x, replace=FALSE)}) 
        if (NSNPS == 1){ hits = t(hits)}
        hdf = gather(data.frame(hits), permute, val)
        print(dim(hdf))
        # Create null catalog:
        prdf = gwdf[hdf$val,]
        prdf$uid = hdf$permute
        prgr = GRanges(paste0('chr', prdf$chrom), IRanges(prdf$chromStart, prdf$chromEnd))
        # Create null intersections:
        qdf = suppressWarnings(data.frame(findOverlaps(prgr, dmgr)))
        qdf$uid = prdf$uid[qdf$queryHits]
        # Non-unique and unique counts:
        qnumdf = aggregate(queryHits ~ uid, qdf, function(x){length(unique(x))})
        qallnumdf = aggregate(queryHits ~ uid, qdf, function(x){length(x)})
        save(prdf, prgr, qdf, qnumdf, qallnumdf, file=gwperfile)
    }
}
