#!/usr/bin/R
# ----------------------------------------------------
# Shuffle full gwas catalogs for significance analysis
# ----------------------------------------------------
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

# -------------------------------------------------
# Only keep GWAS with 10+ lead SNPs (after pruning)
# and with at least 10k individuals.
# -------------------------------------------------
NIND = 10000
NSNP = 10
nsnpdf = aggregate(pValue ~ uid, gwdf, length)
gwssdf = aggregate(sampsize ~ uid, gwssdf, max)
gwssdf = merge(gwssdf, nsnpdf)
keptgw = gwssdf[gwssdf$sampsize >= NIND & gwssdf$pValue >= NSNP,]
kept.uids =  unique(as.character(keptgw$uid))
gwdf = gwdf[gwdf$uid %in% kept.uids,]
TOTSNP = nrow(gwdf)
NCAT = 1000

# Normal intersections:
prgr = GRanges(paste0('chr', gwdf$chrom), IRanges(gwdf$chromStart, gwdf$chromEnd))
gw.qdf = suppressWarnings(data.frame(findOverlaps(prgr, dmgr)))

# Create shuffled catalogs:
for (j in 1:NCAT){
    gwperfile = paste0(perdir, sub(".txt", paste0("_", NCAT, sprintf("_%05d_shuffled_catalog.Rda", j)), gwcatfile))
    if (!file.exists(gwperfile)){
        print(j)
        set.seed(j)
        ind = sample(1:TOTSNP, TOTSNP, replace=FALSE)
        # Create shuffled catalog:
        prdf = gwdf
        prdf$uid = prdf$uid[ind]
        qdf = gw.qdf
        qdf$uid = prdf$uid[qdf$queryHits]
        # Non-unique and unique counts:
        qnumdf = aggregate(queryHits ~ uid, qdf, function(x){length(unique(x))})
        qallnumdf = aggregate(queryHits ~ uid, qdf, function(x){length(x)})
        save(prdf, prgr, qdf, qnumdf, qallnumdf, file=gwperfile)
    }
}
