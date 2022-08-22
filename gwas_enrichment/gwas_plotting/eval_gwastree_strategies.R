#!/usr/bin/R
# ----------------------------------------------------
# Test gwas tree strategies:
# Not parent-child, just weighted averages?
# TODO: Tests:
# - test on good quality (CAD/SCZ/BRCA)
# - Most importantly on bad quality
# - On ones that have the bad nodes (chorion etc.) - flag these...
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
# usetree = 'roadmap'  # For old epi only
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only) }
source(paste0(bindir, 'load_gwastree_analysis.R'))

rm(cdll);
gc()

# ----------------------------------------
# Make sub-directories for data and plots:
# ----------------------------------------
imgdir = paste0(img, "gwas_tree_analysis/test/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'weighted_averages_') 

# Get the leaf reduced information (bag of words, max 3 terms):
leafrep = sapply(1:NN, function(x, max.terms=3){
                     blacklist = c('GLAND', 'TISSUE', 'CELL')
                     x = sub("[0-9]*_","",declist$dec[[x]])
                     # x = unique(x)
                     nx = length(x)
                     if (nx > 1){
                         term_words <- strsplit(x, "[ _,.]");
                         tab_all <- sort(table(unlist(term_words)));
                         tab_all = tab_all[tab_all > 1]
                         tab_all = tab_all[!(names(tab_all) %in% blacklist)]
                         x = paste0(tolower(head(names(sort(tab_all, decreasing=T)), max.terms)), collapse=', ')
                     } else { x = tolower(x) }
                     x = capitalize(x)
                     return(x)})
lind = which(leafrep == '')
nodetissue = nodetissue[order(nodetissue$node),]
leafrep[lind] = nodetissue$GROUP[lind]

inv2x2 = function(A){
    # Determinant
    a <- A[1]
    b <- A[3]
    c <- A[2]
    d <- A[4]
    det <- a * d - b * c
    # Check to see if matrix is singular
    if (det == 0) { stop('Determinant of matrix equals 0, no inverse exists') }
    # Compute inverse matrix elements
    a.inv <- d / det
    b.inv <- -b / det
    c.inv <- -c / det
    d.inv <- a / det
    # Collect the results into a new matrix
    inv.mat <- as.matrix(cbind(c(a.inv,c.inv), c(b.inv,d.inv)))
    return(inv.mat)
}

# --------------------------------------------------------------
# Calculate the aggregate at each node instead of the consensus:
# --------------------------------------------------------------
aggmat.rda.file = 'gwas_tree_analysis/aggregate_bynodes_sparse_matrix.Rda'
if (!file.exists(aggmat.rda.file)){
    # Load enhancer matrix pointers and make sparse matrix:
    mat = read.delim(gzfile(enhmatfile), sep="\t", header=T)
    spmat = sparseMatrix(i=mat$row, j=mat$col)
    # Get the descendants:
    declist$decid = lapply(declist$dec, function(x){labels(dend)[as.numeric(sub("_.*", "", x))]})
    declist$matid = lapply(declist$decid, function(x){which(matnames %in% x)})
    # Transform for aggregation:
    tform = sapply(declist$matid, function(x){(1:833 %in% x)})
    tmat = Matrix(tform * 1)
    aggmat = spmat %*% tmat 
    save(aggmat, file=aggmat.rda.file)
    gc()
} else {
    load(aggmat.rda.file)
}
aggmat = aggmat[enhind,]

# Sweep with matmarg --> do regression
weights = 1 / matmarg[,2]

suid = '29212778 - Coronary artery disease'
# ll = try(get_disjoint_lr(trait=suid, qdf=qdf, type=type, against=against, weights=weights))
# TODO: could model v. marg
# TODO: Try basic reg
# TODO: Try against parent as well.
# TODO: Try both average (agg / nnode) and a weighted average (agg / mmarg  -> frac of enhancer representation)

NE = nrow(enhdf)
pvals = rep(1, NN)
pllk = rep(NA, NN)
dllk = rep(NA, NN)
diffllk = rep(NA, NN)
isnp = rep(0, NN)
# if (return.coeff){ coeffmat = matrix(rep(0, NN * 3), nrow=NN, ncol=3) }

pcoeff = sapply(rep(NA, NN), list)
allsnpset = qdf[qdf$uid == suid, 'subjectHits']
allsnps = rep(0, NE)
allsnps[allsnpset] = 1
NSNP = sum(allsnps)

# Create vectors:
nann = rep(0, NE)
pann = rep(1, NE)
# Simple reg:
pmodel = alt_lr(y=allsnps, x1=pann)
pllk_base = pmodel$loglikelihood

avgmarg = matmarg[,2] / 833
ptm <- proc.time()
pb = txtProgressBar(min=0, max=NN, style = 3)
for (i in 1:NN){
    setTxtProgressBar(pb, i)
    # Sets:
    p = declist$parent[i]
    agg = aggmat[,i]
    # wtd = agg / length(declist$dec[[i]]) # 1. Average
    # wtd = agg * weights # 2. Weighted average
    # wtd = agg / length(declist$dec[[i]]) - avgmarg # 3. Residual of average
    # wtd = agg / length(declist$dec[[i]]) - aggmat[,p] / length(declist$dec[[p]]) # 4. Residual of average vs. parent average
    # wtd = (agg / length(declist$dec[[i]]) - aggmat[,p] / length(declist$dec[[p]])) * weights # 5. Weighted residual of parent 
    wtd = (agg / length(declist$dec[[i]]) - aggmat[,p] / length(declist$dec[[p]])) / avgmarg # 6. Weighted residual of parent 
    # system.time(pmodel <- alt_lr(y=allsnps, x1=wtd))
    # system.time(fastLR(xmat, allsnps))
    # system.time(fastLR(xmat, allsnps, start=c(-6,1)))
    # system.time(fastLR(xmat, allsnps, start=c(-6,0)))
    xmat = cbind(wtd, rep(1, NE))
    out <- fastLR(xmat, allsnps, start=c(.5,-6))
    pcoeff[[i]] = out$coefficients
    dllk[i] = out$loglikelihood
    # Run model (RccpNumerical required):
    # pmodel = alt_lr(y=allsnps, x1=wtd)
    # pcoeff[[i]] = pmodel$coefficients
    # dllk[i] = pmodel$loglikelihood
    diffllk[i] = 2 * (dllk[i] - pllk_base)  # NOTE: the LRT uses 2 x log(R1/R0)
    pvals[i] = pchisq(diffllk[i], df=1, lower.tail=FALSE)
    if (pcoeff[[i]][1] < 0){ pvals[i] = 1 }  # Set 1 if is not positive informative
    cat(i, '\t', round(-log10(pvals[i]),2), '\t', leafrep[i], '\n')
}

# smat = matrix(allsnps)
# xmat = cbind(rep(1, NE), wtd)
## TODO: Make it faster!
## xtx is 1 * 1
#NE 
## 
#xtx = t(xmat) %*% xmat
## Invert:
#xty = t(xmat) %*% smat 
## coeffs = xtx * xty
#xty = aggmat[,i] %*% smat
# Cut off at log10p ~ 20, plot:



# NOTE: Could take all of the matrices, combine.




# TODO: Add stopping + throw NA if looks like it isnt converging.
# calc_disjoint_lr = function(y, x1, x2=NULL, gamma=0.001, max.iter=1500, weights=NULL){
y = allsnps
x1 = wtd
x2 = NULL

NE = length(y)
if (is.null(weights)){ weights = rep(1, NE) } 
if (is.null(x2)){ x2 = rep(0, NE) }
# Count statistics (could be weighted)
N = length(y * weights)
NY = sum(y * weights)
pind = (x1 > 0)
NXP1 = sum((x1 * weights)[pind])
NXN1 = -sum((x1 * weights)[!pind])
NTP1 = sum((x1 * y * weights)[pind])
NFP1 = NXP1 - NTP1 
NFN1 = -sum((x1 * y * weights)[!pind])
NTN1 = NXN1 - NFN1 
NX2 = 0; NTP2 = 0; NFP2 = 0;

# NX2 = sum(x2 * weights)
# NTP2 = sum(x2 * y * weights)
# NFP2 = NX2 - NTP2
# For intercept, not the FP/TP:
NI1 = NY - NTP1 - NTP2 - NFN1 # Y=1, X1=0, X2=0, I=1 (remaining snps)
NI2 = N - NY - NFP1 - NFP2 - NTN1 # Y=0, X1=0, X2=0, I=1 (remaining loc - NY contains NTP1, NTP2)
# Initialize coefficients:
coeff = rep(0,3)
pc = 1e-10 # with pseudo counts
coeff[3] = rsigm((NY - NTP1 - NTP2)/(N - (NX1 + NX2)))
coeff[1] = rsigm((NTP1 + pc)/(NX1)) - coeff[3]
if (NX2 > 0){
    coeff[2] = rsigm((NTP2 + pc)/(NX2)) - coeff[3]
} else {
    coeff[2] = 0
}
j = 1
pctchange = 100
# Update equations::
while (pctchange > 0.0001 && j < max.iter){
    j = j + 1
    oldcoeff = coeff
    # Update x1 (1):
    inner1 = sigm(coeff[1] + coeff[3])
    inner2 = sigm(coeff[2] + coeff[3])
    inner3 = sigm(coeff[3])
    coeff[1] = oldcoeff[1] + 
        gamma * (NTP1 * (1 - inner1) + NFP1 * (0 - inner1))
    # Update x2 (2):
    if (NX2 > 0){
        coeff[2] = oldcoeff[2] +
            gamma * (NTP2 * (1 - inner2) + NFP2 * (0 - inner2))
    }
    # Update intercept (3):
    coeff[3] = oldcoeff[3] + 
        gamma * (NTP1 * (1 - inner1) + NFP1 * (0 - inner1) + 
                 NTP2 * (1 - inner2) + NFP2 * (0 - inner2) + 
                 NI1 * (1 - inner3) + NI2 * (0 - inner3))
    # Percent updated:
    pctchange = sum(abs(coeff - oldcoeff)) / sum(abs(oldcoeff))
}
# Calculate the log likelihood
inner1 = sigm(coeff[1] + coeff[3])
inner2 = sigm(coeff[2] + coeff[3])
inner3 = sigm(coeff[3])
llk = (NTP1 * log(inner1) + NFP1 * log(1 - inner1) + 
       NTP2 * log(inner2) + NFP2 * log(1 - inner2) + 
       NI1 * log(inner3) + NI2 * log(1 - inner3))
# Return the parameters:
return(list(loglikelihood=llk, coefficients=coeff))
# }



# --------------------------------
# Plot the stats as example figure
# --------------------------------

# Under dbdir:
gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
perdir = paste0(gtdir, "permuted_catalogs/")
epref = paste0(usetree, '_e', tol, '_')
regpref = paste0(regdir, epref)
perpref = paste0(perdir, epref)
cmd = paste('mkdir -p', gtdir, regdir, perdir)
system(cmd)

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
suffix = '_adj1000_10.Rda' 
eprefix = paste0(usetree, '_e', tol, '_')

# All Znam in all GWAS past certain 
# for (suid in suidlist) {
# Load regression for suid:
suid = "29212778 - Coronary artery disease" # Have: CAD_UKBIOBANK.gz
plotind = which(uids == suid)
strait = unique(gwdf$trait[gwdf$uid == suid])
spmid = unique(gwdf$pubMedID[gwdf$uid == suid])
sgwssdf = gwssdf[gwssdf$uid == suid,]
ssamp = sgwssdf$sampsize[order(sgwssdf$sampsize, decreasing=T)][1]
sinit = sgwssdf$initSample[order(sgwssdf$sampsize, decreasing=T)][1]
s2init = split.text(sinit, width=90)
ipref = paste0(apref, sprintf("_%05d", plotind))
print(paste(suid, '----- output with prefix', ipref))
regfile = paste0(regpref, ipref, '_lreg', suffix)
load(regfile)

# epigeAt 0.1%
snpfile = paste0(regdir, eprefix, apref, '_epigenomes_hg_all_wsnp', '_adj1000_1.Rda')
load(snpfile)
epi.lp = rmat[suid,]

NTOP=20
# dd$dend = try(set(dd$dend, "branches_lwd", .5))
# Label the top 5 nodes with tissue: 
# ll$rawlp = all.regmat[suid,]
# ll$df$rawlog10p = all.regmat[suid,]
ll$rawlp = -log10(pvals)
ll$df$rawlog10p = -log10(pvals)
ll$log10p = ll$rawlp
ll$log10p[ll$log10p > 12] = 12
ll$df$log10p = ll$log10p 
dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=3, altline=3, altlwd=c(.5, .25), ntop=NTOP)

ldf = ll$df
# ldf = ldf[order(ldf$log10p <= 0, ldf$pout),]
ldf = ldf[order(ldf$log10p <= 0, -ldf$rawlog10p),]
ldf = ldf[1:NTOP,] # TOP N only
MINP = 3
ldf = ldf[ldf$rawlog10p > MINP,]
ntdf = merge(ldf, nodetissue)
ntdf = ntdf[order(ntdf$rawlog10p, decreasing=T),]

ntdf$GROUP = paste0(1:NTOP,'. ', ntdf$GROUP, "\n(", leafrep[ntdf$node], ")")
# ntdf$GROUP = paste0(1:NTOP,'. ', ntdf$GROUP, "\n(", rev(leafrep)[1:NTOP], ")")
ntdf = ntdf[,c('node','GROUP','COLOR')]
names(ntdf) = c('node','symbol','color')
ntdf = merge(ntdf, nodedf)
ntdf$cex=.4

# Setup dendrogram:
# png(paste0(cadpref, 'main_example.png'), res=450, units='in', width=5, height=5)
pdf(paste0('~/test_example.pdf'), width=5, height=5)
# Plot small dendrogram
par(mar=c(0,0,1,0))
set.seed(2) # For legend locations - find a seed that looks ok.
circleplot(dend=dd$dend, lab=lab, fractional=FALSE, hit.track=dd$hitdf, 
           plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=FALSE, 
           nodedf=ntdf, scale=1.5, leaf.pval=epi.lp)
circos.clear()
mtext(strait, side=3, line=0, cex=.7)
mtext(spmid, side=3, line=-.75, cex=.5)
dev.off()

