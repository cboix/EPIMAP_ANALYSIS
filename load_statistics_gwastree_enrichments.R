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
# usetree = 'roadmap'
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions
use.adj = TRUE
use.strict = TRUE
use.onecutoff = FALSE
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
    if (length(statargs) > 6){ use.onecutoff = as.logical(statargs[7]) }
}

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only) }
source(paste0(bindir, 'load_gwastree_analysis.R'))

# if (plotting.only){ rm(dflist) }

# ----------------------------------------
# Make sub-directories for data and plots:
# ----------------------------------------
imgdir = paste0(img, "gwas_tree_analysis/statistics/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
eprefix = paste0(usetree, '_e', tol, '_')
extpref = eprefix
imgpref = paste0(imgdir, eprefix)
treeimgpref = paste0(img, "gwas_tree_analysis/", eprefix)
if (use.adj){ 
    imgpref = paste0(imgpref, 'adj_') 
    treeimgpref = paste0(treeimgpref, 'adj_') 
    extpref = paste0(extpref, 'adj_')
}
if (use.strict){
    imgpref = paste0(imgpref, 'p1_') 
    treeimgpref = paste0(treeimgpref, 'p1_') 
    extpref = paste0(extpref, 'p1_')
}
if (use.onecutoff){
    imgpref = paste0(imgpref, 'onecut_') 
    treeimgpref = paste0(treeimgpref, 'onecut_') 
    extpref = paste0(extpref, 'onecut_')
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

# for (suffix in c("_adj1000_1.Rda", "_adj1000_5.Rda", "_adj1000_10.Rda")){
# Read in regression files:
if (use.adj){ suffix = '_adj1000_10.Rda' } else { suffix = '.Rda' }
if (use.strict){ suffix = '_adj1000_1.Rda' }
if (use.onecutoff){ cutsuff = '_onecut' } else { cutsuff = '' }

all.regfile = paste0(regpref, apref, '_logreg_all', cutsuff, suffix)
print(all.regfile)
if (!file.exists(all.regfile)){
    setpref=''
    # load(paste0(gtdir, 'agg_raw_pvals', setpref, '_20200218.Rda'))
    print("[STATUS] Compiling regression files")
    flist = list.files(path=regdir, pattern=paste0(eprefix, apref, ".*_lreg", suffix))
    uids = sort(as.character(unique(gwdf$uid)))
    all.regmat = matrix(0, nrow=length(uids), ncol=NN)
    rownames(all.regmat) = uids
    for (rfile in flist){
        fnum = sub(paste0("_lreg", suffix), "", sub(paste0(eprefix, apref, "_"), "", rfile))
        fnum = as.numeric(fnum)
        load(paste0(regdir, rfile))
        if (class(ll) == 'list'){
            all.regmat[fnum, ] = -log10(ll$df$pout)
        }
    }
    # Add correction here:
    pcutdf = read.delim(paste0(gtdir, 'cutoffs_l10', setpref, '_all_pernode.tsv'))
    pcutdf = pcutdf[order(pcutdf$node),]
    if (use.adj){
        if (use.onecutoff){
            if (use.strict){
                fper = pcutdf[pcutdf$type == 'All 0.1%','cut']
            } else {
                fper = pcutdf[pcutdf$type == 'All 1%','cut']
            }
            rmatper = all.regmat / fper
        } else {
            if (use.strict){
                fper = pcutdf[pcutdf$type == 'Per 0.1%','cut']
            } else {
                fper = pcutdf[pcutdf$type == 'Per 1%','cut']
            }
            rmatper = sweep(all.regmat, 2, fper, '/')
        }
        all.regmat = 1 * (rmatper > 1) * all.regmat
        # sum(apply(all.regmat > 1, 1, max))
        # sum(apply(all.regmat > 1, 2, max))
    }
    save(all.regmat, file=all.regfile)
} else {
    print("[STATUS] Loading in regression files")
    load(all.regfile)
}
# }

# ------------------------------
# Statistics on the enrichments:
# ------------------------------
cutoff = 3
Z = all.regmat > cutoff
sum(apply(Z, 1, sum) > 0)
sum(apply(Z, 2, sum) > 0)

# Gather into long data.frame:
regwide = data.frame(all.regmat)
colnames(regwide) = 1:ncol(regwide)
regwide$uid = rownames(regwide)
regdf = gather(regwide, node, log10p, -uid)
regdf = regdf[regdf$log10p > cutoff,]
regdf = merge(regdf, data.frame(isleaf=declist$isleaf, 
                                node=1:length(declist$isleaf)))
regdf = merge(regdf, nodetissue)
regdf$type = 'Branch'
regdf$type[regdf$isleaf == 1] = 'Leaf'
gcol = c(colvals$group, 'Multiple' = 'grey80')

# Establish levels/ordering:
odf = rdcol  # Copy don't change rdcol
odf = rbind(odf, c('Multiple','grey80','Other'))
clvs = c('Tissue', 'Primary Cell', 'Immune', 'Neural', 'Stem-like', 'Other')
odf$category = factor(odf$category, levels=clvs)
# Sort alpha then group:
odf = odf[order(odf$GROUP),]
odf = odf[order(odf$category),]
regdf$category = factor(regdf$category, levels=clvs)
regdf$GROUP = factor(regdf$GROUP, levels=odf$GROUP)

# Transform matrices:
tmat = t(sapply(1:NN, function(x){labels(dend3) %in% declist$dec[[x]]}))
ntmat = 1 * t(sapply(1:NN, function(x){odf$GROUP %in% nodetissue$GROUP[nodetissue$node == x]}))
ltmat = 1 * t(sapply(1:NL, function(x){odf$GROUP %in% meta$GROUP[meta$id == labels(dend)[x]]}))
colnames(tmat) = labels(dend3)
colnames(ntmat) = odf$GROUP
colnames(ltmat) = odf$GROUP

# ---------------------------------------------
# Reduce gwas uid to one per trait (587 traits)
# Newest gwas with at least 10k individuals
# And 10+ snps
# ---------------------------------------------
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
keptuids = sapply(traits, function(x){
                      keptgw$uid[keptgw$trait == x][1]})
print(paste0("Kept ", length(keptuids), " traits (newest with sample > ", NIND, ")"))
print(paste0("Kept ", length(keptgwas), " GWAS (newest with sample > ", NIND, ")"))

# NOTE: When loading at some fdr thresh, all below set to 0, 
# so cutoff value could be 0
siguid = names(which(apply(all.regmat > cutoff, 1, sum) > 0))
kuid = siguid[siguid %in% keptuids]
guid = siguid[siguid %in% keptgwas]
print(paste("Kept", length(kuid), "traits (with signif/run)"))
print(paste("Kept", length(guid), "GWAS (with signif/run)"))


# -----------------
# Plot diagnostics:
# -----------------
plot.diagnostics = FALSE
if (plot.diagnostics){
    vimgpref = paste0(img, 'validation/')
    png(paste0(vimgpref, "sampcdf.png"), res=250, width=10, height=8, units='in')
    samp.ord = sort(gwssdf$sampsize)
    samp.ord[samp.ord > 500000] = 500000
    n = length(samp.ord)
    par(mar=c(5,4,1,3))
    plot(samp.ord, (1:n), type = 's', ylim = c(0, n),
         xlab = '', ylab = '', main = '', yaxt='n',xaxt='n')
    cuts = c(1,10,20,30,40,50,100,250,500) * 1000
    par(xpd=FALSE)
    abline(v=cuts, lty='dotted', col=c('black','firebrick'))
    abline(h=hcuts, lty='dotted', col=c('black','firebrick'))
    hcuts = sapply(cuts, function(x){which.min(abs(samp.ord - x))})
    axis(2, at=hcuts, las=1)
    axis(1, at=cuts, las=2)
    dev.off()

    samp.ord = sort(aggregate(sampsize ~ trait, gwssdf, max)$sampsize)
    samp.ord[samp.ord > 500000] = 500000
    png(paste0(vimgpref, "maxtrait_sampcdf.png"), res=250, width=10, height=8, units='in')
    n = length(samp.ord)
    par(mar=c(5,4,1,3))
    plot(samp.ord, (1:n), type = 's', ylim = c(0, n),
         xlab = '', ylab = '', main = '', yaxt='n',xaxt='n')
    cuts = c(1,10,20,30,40,50,100,250,500) * 1000
    par(xpd=FALSE)
    abline(v=cuts, lty='dotted', col=c('black','firebrick'))
    abline(h=hcuts, lty='dotted', col=c('black','firebrick'))
    hcuts = sapply(cuts, function(x){which.min(abs(samp.ord - x))})
    axis(2, at=hcuts, las=1)
    axis(1, at=cuts, las=2)
    dev.off()

    # Match sample size to the rows of regmat:
    ssize = gwssdf$sampsize
    ssize[ssize > 500000] = 500000
    names(ssize) = gwssdf$uid
    ssize = ssize[rownames(all.regmat)]

    # Plot sample size against number of enriched epigenomes at each FDR
    mmarg = apply(all.regmat > 0, 1, sum)
    ind = mmarg > 0
    s1 = ssize[ind]
    s2 = ssize[!ind]
    psuff = paste0(cutsuff, sub("Rda","png", suffix))

    samp.ord = sort(s1)
    samp.ord2 = sort(s2)
    png(paste0(vimgpref, "sampcdf_enrvsnot", psuff), res=250, width=10, height=8, units='in')
    par(mar=c(5,4,1,3))
    n = length(samp.ord)
    plot(samp.ord, (1:n), type = 's', ylim = c(0, max(n,n2)),
         xlab = '', ylab = '', main = '', yaxt='n',xaxt='n')
    n2 = length(samp.ord2)
    lines(samp.ord2, (1:n2), type = 's', ylim = c(0, max(n,n2)), lty='dashed', col='navy')
    cuts = c(1,10,20,30,40,50,100,250,500) * 1000
    par(xpd=FALSE)
    abline(v=cuts, lty='dotted', col=c('black','firebrick'))
    abline(h=hcuts, lty='dotted', col=c('black','firebrick'))
    hcuts = sapply(cuts, function(x){which.min(abs(samp.ord - x))})
    axis(2, at=c(hcuts, max(n,n2)), las=1)
    axis(1, at=cuts, las=2)
    dev.off()

    # Plot number of hits per epigenome vs, its sample size at each FDR
    png(paste0(vimgpref, "sampsize_enrnodes", psuff), res=250, width=10, height=8, units='in')
    par(mar=c(5,4,1,3))
    plot(ssize[ind], mmarg[ind], pch=19, col=rgb(0,0,0,.5), cex=.5, xaxt='n',
         xlab='', ylab='Number of significant nodes (of 1665)')
    axis(1, at=cuts, las=2)
    dev.off()

    # Bin:
    sbins <- cut(ssize, c(0, cuts), include.lowest=T)
    mbins <- cut(mmarg, c(0,1,2,5,10,50,100,200, max(mmarg) + 1), include.lowest=T, right=F)
    mbins = factor(mbins, levels = rev(levels(mbins)))
    bdf = data.frame(size=sbins, n.enr.nodes=mbins)

    gp = ggplot(bdf, aes(size, fill=n.enr.nodes)) + 
        geom_bar(position='fill') + theme_pubclean() + 
        labs(x='Sample size range',y='Fraction') + 
        scale_fill_viridis_d()
    ggsave(paste0('sampsize_binned_enrnodes', psuff), gp, dpi=250, width=9, height=5, units='in')

    # NOTE: There must be an interaction with NSNPs - TODO: plot with snps:
    nsnpdf = aggregate(pValue ~ uid, gwdf, length)
}
