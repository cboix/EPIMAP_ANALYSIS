#!/usr/bin/R
# ---------------------------------------------------
# Count the unique enriched SNPs for all four methods
# ---------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))

# -----------------------
# Defaults for gwas data:
# -----------------------
sets = c('enh', 'rdm', 'mod', 'epi')
use.set = 'epi'

if (use.set == 'rdm'){
    usetree = 'roadmap'
} else {
    usetree = 'enhancers'
}
tol = 2500
singlematch = FALSE
plotting.only = FALSE
runall = TRUE

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only) }
source(paste0(bindir, 'load_gwastree_analysis.R'))
rm(dflist)

# Determine methodology:
use.tree = TRUE
if (use.set == 'epi'){
    flatset = 'epigenomes'
    use.tree = FALSE
} else if (use.set == 'mod'){
    flatset = 'modules'
    use.tree = FALSE
}

# -----------------------------------
# Sub-directories for data and plots:
# -----------------------------------
gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
epref = paste0(usetree, '_e', tol, '_')
regpref = paste0(regdir, epref)
cmd = paste('mkdir -p', gtdir, regdir)
system(cmd)

# ------------------------------------
# Determine which uids we will run on:
# ------------------------------------
uids = sort(as.character(unique(gwdf$uid)))
CHUNKSIZE = length(uids)
# chunk = as.integer(as.numeric(Sys.getenv("SGE_TASK_ID")))
chunk = 1
plotrange = ((chunk - 1)* CHUNKSIZE + 1):min(chunk * CHUNKSIZE, length(uids))
plotuids = uids[plotrange]
# print("[STATUS] Will calculate the following uids:") 
# print(plotuids)

# --------------------------
# Run analysis for each uid:
# --------------------------
# Arguments:
NCAT=1000
type = 'cons'
against = 'parent'
weighted = FALSE
apref = paste0(type, '_', against)
if (weighted){ 
    weights = sqrt(1 / matmarg[,2])
    apref = paste0(apref, '_weighted')
} else {
    weights = NULL
}

# Counts of snps:
nsnpdf = aggregate(pValue ~ uid, gwdf, length)

# Load in the enhancer sets:
if (use.tree){
    enhsets = cdll$diff
    if (use.set == 'roadmap'){
        # Reduce the qdf file to just roadmap enh:
        rdm.enh = enhmap[sort(unique(unlist(enhsets)))]
        qdf = qdf[qdf$subjectHits %in% rdm.enh,]
    }
} else {
    enhsetfile = paste0(flatset, '_enhancer_sets.Rda')
    load(enhsetfile)
    NF = length(enhsets)
}

for (plotind in plotrange){
    cat(plotind, '\n')
    suid = uids[plotind]
    if (use.tree){
        ipref = paste0(apref, sprintf("_%05d", plotind))
    } else {
        ipref = paste0(apref, '_', flatset, sprintf("_%05d", plotind))
    }
    print(paste(plotind, '-', suid, '----- output with prefix', ipref))
    nsnps = nsnpdf$pValue[nsnpdf$uid == suid]
    # Count overlaps at different levels:
    for (lvl in c(1,5,10)){
        if (use.tree){
            cregfile = paste0(regpref, ipref, '_lreg_adj', NCAT, '_', lvl, '.Rda')
        } else {
            cregfile = paste0(regpref, ipref, '_hg_adj', NCAT, '_', lvl, '.Rda')
        }
        if (file.exists(cregfile)){
            load(cregfile)
            if (class(ll) == 'list'){
                # evaluate which are significant (and above cutoff of 2)
                rind = which(ll$rawlp > 2)
                if (length(rind) > 0){
                    nodes = ll$df$node[rind]
                    # Get enhancers in nodes/sets:
                    out = sapply(nodes, function(x){enhsets[[x]]})
                    enhs = sort(enhmap[unique(unlist(out))])
                    # Get overlapping unique and all
                    sub.qdf = qdf[qdf$uid == suid, ]
                    sub.qdf$inenh = 0
                    sub.qdf$inenh[sub.qdf$subjectHits %in% enhs] = 1
                    # Add total unique + all:
                    ll$n.int = sum(sub.qdf$inenh)
                    ll$tot.int = length(sub.qdf$queryHits)
                    ll$n.uniq = length(unique(sub.qdf$queryHits[sub.qdf$inenh == 1]))
                    ll$tot.uniq = length(unique(sub.qdf$queryHits))
                    save(ll, file=cregfile)
                }
            }
        }
    }
}

print(paste("[STATUS] Finished calculating all uids for chunk", chunk))
