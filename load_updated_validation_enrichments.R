#!/usr/bin/R
# -----------------------------------------
# Plot the gwas tree enrichment statistics:
# GWAS-centric analyses/plots:
# -----------------------------------------
statargs=(commandArgs(TRUE))
print(statargs)
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
library(ggplot2)
library(ggpubr)

# Arguments for loading data:
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions
use.adj = TRUE
# use.strict = TRUE
use.strict = FALSE
if (length(statargs)==0) {
    print("Using default arguments. Only loading what is needed for plotting")
} else {        
    usetree = statargs[1]
    tol = as.integer(statargs[2])
    singlematch = as.logical(statargs[3])
    if (length(statargs) > 3){ 
        plotting.only = as.logical(statargs[4])
    }
    if (length(statargs) > 4){ use.adj = as.logical(statargs[5]) }
    if (length(statargs) > 5){ use.strict = as.logical(statargs[6]) }
}

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only) }
source(paste0(bindir, 'load_gwastree_analysis.R'))

rm(dflist)



# -----------------------------------------------
# GWAS with at least 10k individuals and 10+ snps
# -----------------------------------------------
NIND = 10000 # Remove all the low GWAS
NSNP = 10
nsnpdf = aggregate(pValue ~ uid, gwdf, length)
# NOTE: Some duplicates exist.
gwssdf2 = aggregate(sampsize ~ uid, gwssdf, max)
gwssdf = merge(gwssdf, gwssdf2)
gwssdf = unique(merge(gwssdf, nsnpdf))
# 803 GWAS: 
keptgw = gwssdf[gwssdf$sampsize >= NIND & gwssdf$pValue >= NSNP,]
keptgw = keptgw[order(keptgw$pubDate, decreasing=T),]
keptgw = keptgw[order(keptgw$sampsize, decreasing=T),]
traits = unique(keptgw$trait)  # Only unique traits
keptgwas = unique(sort(keptgw$uid)) # All GWAS


# ----------------------------------------
# Make sub-directories for data and plots:
# ----------------------------------------
imgdir = paste0(img, "gwas_tree_analysis/statistics/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
eprefix = paste0(usetree, '_e', tol, '_')
rprefix = paste0('roadmap', '_e', tol, '_')
imgpref = paste0(imgdir, eprefix)
treeimgpref = paste0(img, "gwas_tree_analysis/", eprefix)
if (use.adj){ 
    imgpref = paste0(imgpref, 'adj_') 
    treeimgpref = paste0(treeimgpref, 'adj_') 
}
if (use.strict){
    imgpref = paste0(imgpref, 'p1_') 
    treeimgpref = paste0(treeimgpref, 'p1_') 
}

# Under dbdir:
gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
perdir = paste0(gtdir, "permuted_catalogs/")
epref = paste0(usetree, '_e', tol, '_')
regpref = paste0(regdir, epref)
perpref = paste0(perdir, epref)
cmd = paste('mkdir -p', gtdir, regdir, perdir)
system(cmd)

# ------------------------
# Load in the regressions:
# ------------------------
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

if (use.adj){ suffix = '_adj1000_10.Rda' } else { suffix = '.Rda' }
if (use.strict){ suffix = '_adj1000_1.Rda' }

# Load regression lp and snps for 4 different types of runs:
snpfiles = list()
snpfiles[['enh']] = paste0(regdir, eprefix, apref, '_logreg_all_wsnp', suffix)
snpfiles[['rdm']] = paste0(regdir, rprefix, apref, '_logreg_all_wsnp', suffix)
snpfiles[['mod']] = paste0(regdir, eprefix, apref, '_modules_hg_all_wsnp', suffix)
snpfiles[['epi']] = paste0(regdir, eprefix, apref, '_epigenomes_hg_all_wsnp', suffix)

flist = list()
flist[['enh']] = list.files(path=regdir, pattern=paste0(eprefix, apref, ".*_lreg", suffix))
flist[['rdm']] = list.files(path=regdir, pattern=paste0(rprefix, apref, ".*_lreg", suffix))
flist[['mod']] = list.files(path=regdir, pattern=paste0(eprefix, apref, "_modules.*_hg", suffix))
flist[['epi']] = list.files(path=regdir, pattern=paste0(eprefix, apref, "_epigenomes.*_hg", suffix))

# Cutoffs:
cutdf = read.delim(paste0(gtdir, 'cutoffs_l10_all_pernode.tsv'))
cutdf$set = 'enh'
cdf = read.delim(paste0(gtdir, 'cutoffs_l10_modules_all_pernode.tsv'))
cdf$set = 'mod'
cutdf = rbind(cutdf, cdf)
cutdf = cutdf[cutdf$node == 0,]


# ------------------------------------
# Load raw p-vals, without correction:
# ------------------------------------
snpmats = list()
regmats = list()
nlists = list()
setprefs = c('_', '_modules_','_epigenomes_')
dsetmap = c('enh','mod','epi')
names(dsetmap) = setprefs
for (setpref in setprefs){
    aggpvfile = paste0(gtdir, 'agg_raw_pvals_snps', setpref, '20200218.Rda')
    dset = dsetmap[setpref]
    if (!file.exists(aggpvfile)){
        rmat = matrix(0, nrow=NUID, ncol=nclist[dset])
        smat = matrix(0, nrow=NUID, ncol=nclist[dset])
        rownames(rmat) = uids
        rownames(smat) = uids
        for (rfile in flist[[dset]]){
            if (setpref != '_'){
                fnum = sub(paste0("^.*", apref, ".*_"), "", sub("_hg.*Rda", "", rfile))
            } else {
                fnum = sub(paste0("^.*", apref, ".*_"), "", sub("_lreg.*Rda", "", rfile))
            }
            fnum = as.numeric(fnum)
            load(paste0(regdir, rfile))
            if (class(ll) == 'list'){
                rmat[fnum, ] = -log10(ll$df$pout)
                smat[fnum, ] = ll$isnp
            }
        }

        # Sets of enhancers:
        if (dset == 'enh'){
            enhsets = cdll$diff
        } else if (dset == 'epi'){
            flatset = 'epigenomes'
            enhsetfile = paste0(flatset, '_enhancer_sets.Rda')
            load(enhsetfile)
            NF = length(enhsets)
        } else if (dset == 'mod'){
            flatset = 'modules'
            enhsetfile = paste0(flatset, '_enhancer_sets.Rda')
            load(enhsetfile)
            NF = length(enhsets)
        }
        if (dset == 'epi'){
            load('gwas_hg_stats/ENH_wH3K27ac100_hg_processed_enr.Rda')
        }

        rmatlist = list()
        smatlist = list()
        nhitlist = list()
        for (lvl in c('1%','0.1%')){
            # Reduce to kept only:
            if (dset == 'epi'){
                if (lvl == '1%'){
                    fdr = '99%' 
                } else { 
                    fdr = '99.9%' 
                }
                regmat = calist[[fdr]][keptgwas,]
                snpmat = smat[keptgwas,]
                snpmat = snpmat * (regmat > 0)
            } else {
                cutoff = cutdf[cutdf$set == dset & cutdf$type == paste('All',lvl), 'cut']
                regmat = rmat * (rmat > cutoff)
                snpmat = smat * (rmat > cutoff)
            }

            regmat = regmat[keptgwas,]
            snpmat = snpmat[keptgwas,]
            sum(apply(regmat, 1, max) > 0) # Number passing
            nlist = matrix(0, nrow(regmat), 4)
            rownames(nlist) = rownames(regmat)
            colnames(nlist) = c('n.int','tot.int','n.uniq','tot.uniq')

            # Count the SNPs in these:
            for (i in 1:nrow(regmat)){
                cat(i, '\n')
                suid = rownames(regmat)[i]
                rind = which(regmat[suid,] > 0)
                if (length(rind) > 0){
                    nodes = ll$df$node[rind]
                    # Get enhancers in nodes/sets:
                    out = sapply(rind, function(x){enhsets[[x]]})
                    enhs = sort(enhmap[unique(unlist(out))])
                    # Get overlapping unique and all
                    sub.qdf = qdf[qdf$uid == suid, ]
                    sub.qdf$inenh = 0
                    sub.qdf$inenh[sub.qdf$subjectHits %in% enhs] = 1
                    # Add total unique + all:
                    nlist[suid, 'n.int'] = sum(sub.qdf$inenh)
                    nlist[suid, 'tot.int'] = length(sub.qdf$queryHits)
                    nlist[suid, 'n.uniq'] = length(unique(sub.qdf$queryHits[sub.qdf$inenh == 1]))
                    nlist[suid, 'tot.uniq']= length(unique(sub.qdf$queryHits))
                }
            }
            print(colSums(nlist))
            rmatlist[[lvl]] = regmat
            smatlist[[lvl]] = snpmat
            nhitlist[[lvl]] = nlist
        }
        save(rmatlist, smatlist, nhitlist, file=aggpvfile)
    } else {
        load(aggpvfile)
    }
    # Add to the full set:


}

sets = c('enh','rdm','mod','epi')
nlist = NULL
regmats = list()
snpmats = list()
nlists = list()
nclist = c(NN, 417, 300, NL)
names(nclist) = sets
runall = FALSE
for (dset in sets){
    if (!file.exists(snpfiles[[dset]]) || runall){
        print(paste("[STATUS] Compiling regression files -", dset))
        uids = sort(as.character(unique(gwdf$uid)))
        NUID = length(uids)
        rmat = matrix(0, nrow=NUID, ncol=nclist[dset])
        smat = matrix(0, nrow=NUID, ncol=nclist[dset])
        nlist = list(n.int = rep(0, NUID), tot.int = rep(0, NUID),
                     n.uniq = rep(0, NUID), tot.uniq = rep(0, NUID))
        rownames(rmat) = uids
        rownames(smat) = uids
        for (rfile in flist[[dset]]){
            if (dset %in% c('epi','mod')){
                fnum = sub(paste0("^.*", apref, ".*_"), "", sub("_hg.*Rda", "", rfile))
            } else {
                fnum = sub(paste0("^.*", apref, ".*_"), "", sub("_lreg.*Rda", "", rfile))
            }
            fnum = as.numeric(fnum)
            load(paste0(regdir, rfile))
            if (class(ll) == 'list'){
                rmat[fnum, ] = ll$rawlp
                smat[fnum, ] = ll$isnp
                if (!is.null(ll$n.int)){
                    nlist$n.int[fnum] = ll$n.int 
                    nlist$tot.int[fnum] = ll$tot.int 
                    nlist$n.uniq[fnum] = ll$n.uniq
                    nlist$tot.uniq[fnum] = ll$tot.uniq
                }
            }
        }
        rmat[is.na(rmat)] = 0
        save(rmat, smat, nlist, file=snpfiles[[dset]])
    } else {
        print("[STATUS] Loading in regression files")
        load(snpfiles[[dset]])
    }
    regmats[[dset]] = rmat
    snpmats[[dset]] = smat
    nlists[[dset]] = nlist
}


NIND = 20000 # Remove all the really low GWAS
keptgw = gwssdf[gwssdf$sampsize > NIND,]

for (dset in sets){
    rmat = regmats[[dset]] > 0
    rmat[is.na(rmat)] = 0
    kuid = rownames(rmat)
    kuid = kuid[kuid %in% keptgw$uid]
    print(length(kuid))
    print(paste(dset, sum(apply(rmat, 1, sum) > 0)))
    print(paste(dset, sum(apply(rmat[kuid,], 1, sum) > 0)))
}

# Look at co-top sets:
rmat = regmats[['epi']]
rmat = regmats[['mod']]
rmat[is.na(rmat)] = 0
sum(rmat > 0)

kuid = names(which(apply(rmat > 0, 1, sum) > 0))
sum(apply(rmat > 0, 1, sum) > 0)
sum(apply(rmat > 0, 2, sum) > 0)
kuid = kuid[kuid %in% keptgw$uid]

# Look at group sets:
kmat = rmat[kuid,]
kassign = apply(kmat, 1, which.max)
lk = lapply((1:833)-1, function(x){names(which(kassign == x))})
lind = which(lapply(lk, length) > 1)
# print(lind)
# for (li in lind){
#     print("---------")
#     print(lk[[li]])
# }
# apply(rmat[kuid,] > 0, 2, sum)
# sort(apply(rmat[kuid,] > 10, 2, sum))

# -------------------------
# Collect test set lengths:
# -------------------------
lensfile = 'consensus_object_lengths_all_main_types.Rdata'
if (!file.exists(lensfile)){
    lens = list()
    cdlenfile = paste0('consensus_object_lengths_', 'enhancers', '_062819.Rdata')
    load(cdlenfile)
    lens[['enh']] = cdlenlist$diff
    cdlenfile = paste0('consensus_object_lengths_', 'roadmap', '_062819.Rdata')
    load(cdlenfile)
    lens[['rdm']] = cdlenlist$diff
    enhsetfile = paste0('modules_enhancer_sets.Rda')
    load(enhsetfile)
    lens[['mod']] = lenlist
    enhsetfile = paste0('epigenomes_enhancer_sets.Rda')
    load(enhsetfile)
    lens[['epi']] = lenlist
    save(lens, file=lensfile)
} else {
    load(lensfile)
}

