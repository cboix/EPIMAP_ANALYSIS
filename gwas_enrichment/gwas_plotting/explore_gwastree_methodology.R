#!/usr/bin/R
# -------------------------------------
# Exploration + verification of results
# under different methodologies
# -------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
library(glmnet)
library(igraph)

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

rm(dflist)

# ----------------------------------------
# Make sub-directories for data and plots:
# ----------------------------------------
imgdir = paste0(img, "gwas_tree_analysis/examples/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
treeimgpref = paste0(imgdir, usetree, '_e', tol, '_')

# Under dbdir:
gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
perdir = paste0(gtdir, "permuted_catalogs/")
epref = paste0(usetree, '_e', tol, '_')
regpref = paste0(regdir, epref)
perpref = paste0(perdir, epref)
cmd = paste('mkdir -p', gtdir, regdir, perdir)
system(cmd)

# Preparation:
tree.enhsets = cdll$diff
tree.enhlens = cdlenlist$diff
tree.par.enhsets = cdll$cons
tree.par.enhlens = cdlenlist$cons
nodetissue = nodetissue[order(nodetissue$node),]

# Plot the independent log pvalues against the tree ones:
plot.small.tree = function(ll, suid, minp=3, title=''){
    plotind = which(uids == suid)
    strait = unique(gwdf$trait[gwdf$uid == suid])
    spmid = unique(gwdf$pubMedID[gwdf$uid == suid])
    sgwssdf = gwssdf[gwssdf$uid == suid,]
    sinit = sgwssdf$initSample[order(sgwssdf$sampsize, decreasing=T)][1]
    s2init = split.text(sinit, width=90)
    # Setup dendrogram:
    dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=minp, altline=0)
    dd$dend = try(set(dd$dend, "branches_lwd", .5))
    # Label the top 5 nodes with tissue: 
    NTOP=5
    ldf = ll$df
    ldf = ldf[order(ldf$log10p <= 0, ldf$pout),]
    ldf = ldf[1:NTOP,] # TOP N only
    ldf = ldf[ldf$rawlog10p > minp,]
    ntdf = merge(ldf, nodetissue)
    ntdf = ntdf[,c('node','GROUP','COLOR')]
    names(ntdf) = c('node','symbol','color')
    ntdf = merge(ntdf, nodedf)
    ntdf$cex=.5
    # Plot small dendrogram
    par(mar=c(1,0,1,0))
    circleplot(dend=dd$dend, lab=lab, fractional=FALSE, hit.track=dd$hitdf, 
               plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=TRUE, nodedf=ntdf, scale=2)
    circos.clear()
    if (title == ''){
        title = strait 
    } else {
        title = paste(title, " - ", strait)
    }
    mtext(title, side=3, line=0, cex=.7)
    mtext(spmid, side=3, line=-.75, cex=.5)
    mtext(s2init, side=1, line=0, cex=.3)
}


plot.compare = function(x, y, xlab, ylab, main=''){
    min.lp = apply(cbind(x, y), 1, min)
    max.lp = apply(cbind(x, y), 1, max)
    lims = c(min(min.lp), max(max.lp))
    plot(x, y, pch=21, 
         cex = apply(cbind(x/10, rep(.25, NN)), 1, max), 
         bg = sapply(nodetissue$COLOR, tsp.col), lwd=min.lp > 10,
         ylim=lims, xlim=lims, 
         xlab=xlab, ylab=ylab, main=main)
    ngroup = nodetissue$GROUP
    ngroup[x < 10] = ''
    text(x, y, ngroup, col=nodetissue$COLOR, cex=x/100)
    abline(h=0)
    abline(v=0)
    abline(0,1, lty=2)
}

# ----------------------------------
# Determine which uids we will plot:
# ----------------------------------
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
uids = sort(as.character(unique(gwdf$uid)))
suid = "29212778 - Coronary artery disease"
plotind = which(uids == suid)
sub.qdf = qdf[qdf$uid == suid, ]

# Load current with correction:
ipref = paste0(apref, sprintf("_%05d", plotind))
regfile = paste0(regpref, ipref, '_lreg.Rda')
suffix = '_adj1000_1.Rda'
cregfile = paste0(regpref, ipref, '_lreg', suffix)
load(regfile)
fixed.ll = ll
load(cregfile)
corr.ll = ll

# ---------------------
# ORIGINAL METHODOLOGY:
# Using verbose version (return coeff).
# ---------------------
recalc.ll = get_disjoint_lr(trait=suid, qdf=qdf, type=type, against=against, 
                            weights=weights, verbose=T, return.coeff=T, 
                            max.iter=1500)

cmat = recalc.ll$coeff
lims = c(-1, max(cmat[,c(1:2)]))
png('~/test_coeff.png', res=450, units='in', width=9, height=9)
lp.par = recalc.ll$rawlp[declist$parent]
plot(cmat[,2], cmat[,1], pch=21, 
     cex = apply(cbind(recalc.ll$rawlp /10, rep(.25, NN)), 1, max), 
     bg = sapply(nodetissue$COLOR, tsp.col), lwd=lp.par > 10,
     ylim=lims, xlim=lims,
     ylab='Parent Coefficient',
     xlab='Node (diff) Coefficient')
abline(h=0)
abline(v=0)
abline(0,1, lty=2)
dev.off()

png('~/test_fixedpval.png', res=450, units='in', width=6, height=6)
plot(fixed.ll$rawlp, recalc.ll$rawlp, pch=19, 
     xlab='Original', ylab='Re-calculated')
abline(0,1)
dev.off()

if (is.null(verif.ll)){
    verif.ll = get_disjoint_lr(trait=suid, qdf=qdf, type=type, against=against, 
                               verbose=T, return.coeff=T, alt.lr=TRUE)
}


png('~/comp_verification.png', res=450, units='in', width=6, height=6)
plot(recalc.ll$rawlp, verif.ll$rawlp, pch=19, 
     xlab='Re-calculated', ylab='Verification - External Method')
abline(0,1)
dev.off()

union.ll = get_disjoint_lr(trait=suid, qdf=qdf, type='novel', against='parent', 
                            weights=weights, verbose=T, return.coeff=T)

cmat = union.ll$coeff
lims = c(-1, max(cmat[,c(1:2)]))
lims = range(cmat[,c(1:2)])
png('~/test_coeff_union.png', res=450, units='in', width=9, height=9)
cmat = union.ll$coeff
lp.par = union.ll$rawlp[declist$parent]
plot(cmat[,2], cmat[,1], pch=21, 
     cex = apply(cbind(union.ll$rawlp /10, rep(.25, NN)), 1, max), 
     bg = sapply(nodetissue$COLOR, tsp.col), lwd=lp.par > 10,
     ylim=lims, xlim=lims,
     ylab='Parent Coefficient',
     xlab='Node (diff) Coefficient')
ngroup = nodetissue$GROUP
ngroup[union.ll$rawlp < 10] = ''
text(cmat[,2], cmat[,1], ngroup, 
     col=nodetissue$COLOR, cex=union.ll$rawlp/100)
abline(h=0)
abline(v=0)
abline(0,1, lty=2)
dev.off()

diff.ll = get_disjoint_lr(trait=suid, qdf=qdf, type='diff', against='parent', 
                            weights=weights, verbose=T, return.coeff=T)

cmat = diff.ll$coeff
lims = range(cmat[,c(1:2)])
lims = c(-1, max(cmat[,c(1:2)]))
png('~/test_coeff_diff.png', res=450, units='in', width=9, height=9)
cmat = diff.ll$coeff
lp.par = diff.ll$rawlp[declist$parent]
plot(cmat[,2], cmat[,1], pch=21, 
     cex = apply(cbind(diff.ll$rawlp /10, rep(.25, NN)), 1, max), 
     bg = sapply(nodetissue$COLOR, tsp.col), lwd=lp.par > 10,
     ylim=lims, xlim=lims,
     ylab='Parent Coefficient',
     xlab='Node (diff) Coefficient')
ngroup = nodetissue$GROUP
ngroup[diff.ll$rawlp < 10] = ''
text(cmat[,2], cmat[,1], ngroup, 
     col=nodetissue$COLOR, cex=diff.ll$rawlp/100)
abline(h=0)
abline(v=0)
abline(0,1, lty=2)
dev.off()


# ----------------------------
# Also evaluate independently:
# ----------------------------
diff.dflist = eval.intersections(cdll$diff, sub.qdf)
novel.dflist = eval.intersections(cdll$novel, sub.qdf)
idiff.ll = get_pvalues_flat(suid, diff.dflist, cdlenlist$diff, sub.qdf, NF=NN, cutp=12)
inovel.ill = get_pvalues_flat(suid, novel.dflist, cdlenlist$novel, sub.qdf, NF=NN, cutp=12)

png('~/comp_cons_recalc_indpt.png', res=450, units='in', width=9, height=9)
plot.compare(x=recalc.ll$rawlp, y=idiff.ll$rawlp, 
             xlab='Tree Regression log10p',
             ylab='Independent log10p (hyper geom)',
             main='Consensus difference - method comparsion\n(black border if both signif.)')
dev.off()

png('~/comp_union_recalc_indpt.png', res=450, units='in', width=9, height=9)
plot.compare(x=union.ll$rawlp, y=inovel.ill$rawlp, 
             xlab='Tree Regression log10p',
             ylab='Independent log10p (hyper geom)',
             main='Union difference - method comparsion\n(black border if both signif.)')
dev.off()

png('~/comp_union_cons.png', res=450, units='in', width=9, height=9)
plot.compare(x=recalc.ll$rawlp, y=union.ll$rawlp,
             xlab='Consensus difference log10p',
             ylab='Union difference log10p',
             main='Method comparsion - tree methods (cons vs. union)\n(black border if both signif.)')
dev.off()

png('~/comp_union_cons_indpt.png', res=450, units='in', width=9, height=9)
plot.compare(x=idiff.ll$rawlp, y=inovel.ill$rawlp,
             xlab='Consensus difference log10p',
             ylab='Union difference log10p',
             main='Method comparsion - independent (cons vs. union)\n(black border if both signif.)')
dev.off()


TOPN = 25
png(paste0('~/tree_top', TOPN, '_cons.png'), res=450, units='in', width=4, height=4)
use.ll = recalc.ll
ind = order(use.ll$rawlp, decreasing=T)
# Make sure coeff not negative:
ind = ind[use.ll$coeff[ind, 2] > 0]
ind = ind[1:TOPN]
use.ll$rawlp[-ind] = 0
use.ll$log10p[-ind] = 0
plot.small.tree(use.ll, suid, title='Consensus')
dev.off()


png(paste0('~/tree_top', TOPN, '_union.png'), res=450, units='in', width=4, height=4)
use.ll = union.ll
ind = order(use.ll$rawlp, decreasing=T)
# Make sure coeff not negative:
ind = ind[use.ll$coeff[ind, 2] > 0]
ind = ind[1:TOPN]
use.ll$rawlp[-ind] = 0
use.ll$log10p[-ind] = 0
plot.small.tree(use.ll, suid, title='Union')
dev.off()


png(paste0('~/tree_top', TOPN, '_int.png'), res=450, units='in', width=4, height=4)
min.lp = apply(cbind(recalc.ll$rawlp, union.ll$rawlp), 1, min)
min.lp[min.lp < 3] = 0
ind = order(min.lp, decreasing=T)
use.ll = union.ll
ind = ind[use.ll$coeff[ind, 2] > 0]
ind = ind[1:TOPN]
min.lp[-ind] = 0
use.ll$rawlp = min.lp
use.ll$df$rawlog10p= min.lp
min.lp[min.lp > 12] = 12
use.ll$log10p = min.lp
use.ll$df$log10p = min.lp
plot.small.tree(use.ll, suid, title='Intersection')
dev.off()



use.ll = recalc.ll
plot(corr.ll$rawlp, fixed.ll$rawlp)

# Get the snp sets for each node:
rind = which(use.ll$rawlp > 3)
nodes = use.ll$df$node[rind]
lp = use.ll$rawlp[rind]
parents = declist$parent[nodes]
siblings = sapply(nodes, function(x){
                      p = declist$parent[x]
                      d = which(declist$parent == p)
                      d[d != x] })
lp.par = - log10(use.ll$df$pout)[parents]
pn = parents[parents %in% nodes]
nodegroup = nodetissue$GROUP[nodes]
nodecol = nodetissue$COLOR[nodes]
nodelp = round(use.ll$rawlp[rind],1)
nodenames = paste0(nodegroup, ' - ', nodelp)
names(nodecol) = nodenames

# Get enhancers in nodes/sets:
out = sapply(nodes, function(x){tree.enhsets[[x]]})
omat = sapply(out, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 
out2 = sapply(parents, function(x){tree.enhsets[[x]]})
omat2 = sapply(out2, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 
out3 = sapply(parents, function(x){tree.par.enhsets[[x]]})
omat3 = sapply(out3, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 
out4 = sapply(siblings, function(x){tree.enhsets[[x]]})
omat4 = sapply(out4, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 

r1 = apply(omat,2,sum) / tree.enhlens[nodes]
r2 =  apply(omat2, 2, sum) / tree.enhlens[parents]
r3 =  apply(omat3, 2, sum) / tree.par.enhlens[parents]
r4 =  apply(omat4, 2, sum) / tree.enhlens[siblings]

# png('~/test2.png', res=450, units='in', width=9, height=9)
png('~/test3.png', res=450, units='in', width=9, height=9)
plot(r1, r3, ylim =c(0,max(r1,r3)), xlim=c(0,max(r1,r3)), 
     pch=21, cex = lp/10, bg=sapply(nodecol, tsp.col), lwd = lp.par > 10,
     xlab='# hits in node diff / # enhancers',
     ylab='# hits in parent consensus/ # enhancers', 
     main='Signif. nodes for Coronary Artery Disease\n(size=node log10p, text=parent log10p)')
# text(r1, r3, nodenames, col=nodecol, cex=lp/50)
text(r1, r3, round(lp.par,1), col=nodecol, cex=lp/50)
abline(0,1)
dev.off()



# -----------------------------
# Compare against other traits:
# -----------------------------
plot.scatter.traits = function(df){
    sp=.1
    par(mar=c(2,2, sp,sp))
    xymat = cbind(df$main, df$against)
    min.lp = apply(xymat, 1, min)
    lims = range(c(df$main, df$against))
    plot(df$against, df$main, pch=21, axes=F,
         cex = apply(cbind(df$main/20, rep(.25, NN)), 1, max), 
         bg = sapply(df$COLOR, tsp.col), 
         lwd=.5 * (min.lp > 10), ylim=lims, xlim=lims)
    box()
    if (nchar(as.character(df$trait[1])) > 60){
        mtext(split.text(df$trait[1],20),1,line=.1, cex=.9)
    } else {
        mtext(df$trait[1],1, cex=.9, line=.1)
    }
    mtext(df$main.trait[1], 2, cex=.9, line=.1)
    ngroup = df$GROUP
    ngroup[min.lp < 10] = ''
    text(df$against, df$main, ngroup, col=nodetissue$COLOR, cex=df$main/100)
    abline(0,1, lty=2, lwd=.5)
    # Add regression line + r2.
    model = lm(main ~ against, df)
    abline(model$coeff[1], model$coeff[2], lty='dotted', col='firebrick')
    ar2 = round(summary(model)$adj.r.squared, 3)
    text(x=lims[1] + .99 * diff(lims), 
         y=lims[1] + .01 * diff(lims), adj=1, col='firebrick',
         as.expression(c(bquote(paste("adj. ", R^2 == .(ar2)))))) 
}

suid = "29212778 - Coronary artery disease"
slist = c("29507422 - High density lipoprotein cholesterol levels",
          "29507422 - Low density lipoprotein cholesterol levels",
          "29507422 - Total cholesterol levels",
          "29507422 - Triglycerides")

slist2 = c("26343387 - Myocardial infarction",
           "29403010 - Mean arterial pressure",
           "30595370 - Systolic blood pressure",
           "27841878 - Systolic blood pressure",
           "30595370 - Cardiovascular disease",
           "29892015 - Atrial fibrillation",
           "30061737 - Atrial fibrillation",
           "24952745 - QT interval")


slist3 = c('23969696 - Fibrinogen', 
           '30367059 - Thyroid stimulating hormone',
           '29422604 - Pancreatic cancer',
           '30529582 - Colorectal cancer', 
           # '26831199 - Chronic kidney disease', 
           '30504769 - Gallstone disease', 
           "29531354 - Ischemic stroke",
           # "29531354 - Stroke",
           "27588450 - Glomerular filtration rate",
           # "25628336 - Motion sickness", 
           "30038396 - Cognitive performance (MTAG)",
           "30038396 - Educational attainment (MTAG)",
           "30289880 - Lipid traits (pleiotropy) (HIPO component 1)",
           '23326517 - Age-related macular degeneration',
           "30535121 - Macular thickness",
           "30054594 - Intraocular pressure",
           "26192919 - Ulcerative colitis",
           "27182965 - Allergy",
           "30617256 - Alzheimer's disease or family history of Alzheimer's disease")

odf = rdcol  # Copy don't change rdcol
odf = rbind(odf, c('None','grey80','Other'))
clvs = c('Tissue', 'Primary Cell', 'Immune', 'Neural', 'Stem-like', 'Other')
odf$category = factor(odf$category, levels=clvs)
# Sort alpha then group:
odf = odf[order(odf$GROUP),]
odf = odf[order(odf$category),]

# Read in regression files (for adjusted log10p):
plotind = which(uids == suid)
ipref = paste0(apref, sprintf("_%05d", plotind))
regfile = paste0(regpref, ipref, '_lreg.Rda')
load(regfile)
fixed.ll = ll
lp = fixed.ll$rawlp
# Alternatively, load current with correction?
# suffix = '_adj1000_1.Rda'
# cregfile = paste0(regpref, ipref, '_lreg', suffix)
# load(cregfile)
# corr.ll = ll

vdf = c()
for (suid2 in c(slist, slist2, slist3)){
    # Read in regression files (for adjusted log10p):
    plotind2 = which(uids == suid2)
    ipref = paste0(apref, sprintf("_%05d", plotind2))
    regfile = paste0(regpref, ipref, '_lreg.Rda')
    load(regfile)
    vdf = rbind(vdf, data.frame(main=lp, against=ll$rawlp, node=1:NN, uid=suid2, 
                                main.trait=sub("^[0-9]* - ", "", suid), trait=sub("^[0-9]* - ", "", suid2)))
}
vdf = merge(vdf, nodetissue)



png('~/CAD_vs_lipids.png', res=450, units='in', width=12, height=3)
layout(matrix(c(1:4), 1,4), TRUE)
for (suid2 in slist){
    df = vdf[vdf$uid == suid2,]
    plot.scatter.traits(df)
}
dev.off()

png('~/CAD_vs_others.png', res=450, units='in', width=12, height=6)
layout(matrix(c(1:8), 2,4), TRUE)
for (suid2 in slist2){
    df = vdf[vdf$uid == suid2,]
    plot.scatter.traits(df)
}
dev.off()

png('~/CAD_vs_all.png', res=450, units='in', width=12, height=9)
layout(matrix(c(1:12), 3,4, byrow=T), TRUE)
for (suid2 in c(slist,slist2)){
    df = vdf[vdf$uid == suid2,]
    plot.scatter.traits(df)
}
dev.off()


png('~/CAD_vs_external.png', res=450, units='in', width=12, height=12)
layout(matrix(c(1:16), 4,4, byrow=T), TRUE)
for (suid2 in slist3){
    df = vdf[vdf$uid == suid2,]
    plot.scatter.traits(df)
}
dev.off()


# -------------------------------------------
# Now run regressions on all of the kept GWAS 
# and look at which associate
# more or less with CAD.
# -------------------------------------------
strait = sub("^[0-9]* - ", "", suid)
NIND = 20000 # Remove all the really low GWAS
keptgw = gwssdf[gwssdf$sampsize > NIND,]
keptgwas = unique(keptgw$uid)
NGWAS = length(keptgwas)
coeffmat = matrix(rep(0, 2 * NGWAS), nrow=NGWAS, ncol=2)
lpmat = matrix(rep(0, NN * NGWAS), nrow=NGWAS, ncol=NN)
r2vec = rep(0, NGWAS)
rownames(coeffmat) = keptgwas
rownames(lpmat) = keptgwas
names(r2vec) = keptgwas
coeffs = c()
suffix = "_adj1000_10.Rda"
for (suid2 in keptgwas){
    plotind2 = which(uids == suid2)
    ipref = paste0(apref, sprintf("_%05d", plotind2))
    regfile = paste0(regpref, ipref, '_lreg', suffix)
    load(regfile)
    if (ll != 'error'){
        # 
        lpmat[suid2,] = ll$rawlp
        df = data.frame(main=lp, against=ll$rawlp)
        model = lm(main ~ against, df)
        coeffmat[suid2,] = model$coeff
        r2vec[suid2] = summary(model)$adj.r.squared
    }
}

# NOTE: Load the corrected version!
snpfile = paste0(regdir, epref, apref, '_logreg_all_wsnp', suffix)
load(snpfile)
lpmat = rmat[keptgwas,]

lpmat2 = lpmat
lpmat2[lpmat2 > 12] = 12
dt = dist(lpmat2, 'cosine')
ht = hclust(dt, 'ward.D')
cocl = order.optimal(dt, ht$merge)$order
rn <- names(cocl)[cocl]
dt = dist(t(lpmat2), 'cosine')
ht = hclust(dt, 'ward.D')
cocl = order.optimal(dt, ht$merge)$order

png('~/reord_lpmat.png', res=450, units='in', width=9, height=9)
par(mar=rep(0.1, 4))
image(lpmat2[rn,cocl], axes=F)
dev.off()
# sort(r2vec)

# ----------------------------------
# Run lasso to predict CAD 
# based on the other regressions -
# what component traits sum to CAD?
# ----------------------------------
model = glmnet(t(lpmat[rownames(lpmat) != suid,]), lp, alpha=.8)
summary(model)
plot(model, xvar = "dev", label = TRUE)
cf = as.matrix(coefficients(model))

ord = order(apply(cf > 0, 1, sum))
image(cf[ord,] > 0)
rownames(cf)[ord]


# -------------------------------------------------
# As an experiment, run all pairs of regressions
# NOTE: will need to correct for promiscuous traits.
# -------------------------------------------------
gather.sqmat = function(mat, names){
    rn = rownames(mat)
    mat = data.frame(mat)
    colnames(mat) = rn
    rownames(mat) = rn
    mat$G1 = rn
    df = gather(mat, G2, metric, -G1)
    names(df) = names
    return(df)
}

# Remove zero-regs (not greater than 3, adjusted):
lind = which(apply(lpmat > 3, 1 ,sum) > 0)
red.lpmat = lpmat[lind,]
compfile = paste0(epref, 'comp_regs_padjust', suffix)
slist = rownames(red.lpmat)
if(!file.exists(compfile)){
    spearman = cor(t(red.lpmat), method='spearman')
    pearson = cor(t(red.lpmat), method='pearson')
    rsqmat = pearson
    rsqmat[] = 0
    i = 0
    for (suid1 in slist){
        i = i + 1
        print(i)
        lp1 = red.lpmat[suid1,]
        for (suid2 in slist){
            df = data.frame(main=red.lpmat[suid1,], against=red.lpmat[suid2,])
            model = lm(main ~ against, df)
            rsqmat[suid1, suid2] = summary(model)$adj.r.squared
        }
    }
    # Turn into compdf:
    rdf = gather.sqmat(rsqmat, c('G1','G2','rsq'))
    pdf = gather.sqmat(pearson, c('G1','G2','pearson'))
    sdf = gather.sqmat(spearman, c('G1','G2','spearman'))
    compdf = merge(merge(rdf, pdf), sdf)
    save(compdf, rsqmat, pearson, spearman, file=compfile)
} else {
    load(compfile)
}


# --------------------------------------
# Preliminary exploration of statistics:
# - are there any negative correlations?
# - what are the strongest correlations?
# - plots of different measure matrices.
# --------------------------------------
summary(compdf)
# tail(compdf[compdf$spearman < -0.1,])


# --------------------------------------
# Plot a network for a specific measure:
# --------------------------------------
# TODO: Do we need to calibrate the measure? 
measure = 'rsq'
cutoff = .1 # Magnitude cutoff for rsq or other measure
zcutoff = 2 # Magnitude cutoff for rsq or other measure
by.zscore=TRUE
ddf = compdf[,c('G1','G2',measure)]
names(ddf)[3] = 'measure'
ddf[ddf$G1 != ddf$G2,]
mdf = aggregate(measure ~ G1, ddf, mean)
sdf = aggregate(measure ~ G1, ddf, sd)
names(mdf)[2] = 'mean'
names(sdf)[2] = 'sd'
ddf = merge(merge(ddf, mdf), sdf)
ddf$zscore = (ddf$measure - ddf$mean) / ddf$sd
ddf$is.neg = ddf$measure < 0  # Express negative edges as dashed lines.
if(by.zscore){
    ddf = ddf[ddf$zscore >= zcutoff,]
} else {
    ddf = ddf[abs(ddf$measure) >= cutoff,]
}
# TODO: Find way to reduce edges in rsq formulation.
# ddf = ddf[ddf$G1 != ddf$G2,]
nodes = sort(unique(c(as.character(ddf$G1), 
                      as.character(ddf$G2))))

# Keep edges in opposite direction (diff meaning for rsq)
ddf$G1 = factor(ddf$G1, levels=nodes)
ddf$G2 = factor(ddf$G2, levels=nodes)
if (measure != 'rsq'){
    ddf = ddf[as.numeric(ddf$G1) < as.numeric(ddf$G2),]
}
ddf = ddf[as.numeric(ddf$G1) < as.numeric(ddf$G2),]
print(paste0("Kept ", dim(ddf)[1], " edges of ", 
             dim(compdf)[1], " total edges (", 
             round(dim(ddf)[1] / dim(compdf)[1] * 100, 2), "%)"))

# Get maximal node:
tform = make.tform(nodetissue$GROUP, u=odf$GROUP)
colnames(tform) = odf$GROUP
slist = rownames(red.lpmat)
uslist = unique(slist) 
lpmat.group = red.lpmat[nodes,] %*% tform 
m = apply(ddf, 1, function(x){which.max(lpmat.group[x[1],] * lpmat.group[x[2],])})
ddf$GROUP = colnames(lpmat.group)[m]
ddf = merge(ddf, odf)
ddf = ddf[,c('G1','G2','GROUP','measure','is.neg','COLOR','category', 'zscore')]

by.max=FALSE
if (by.max){
    diagcto = apply(red.lpmat, 1, which.max)
    ctodf = data.frame(gw.node=names(diagcto), node=diagcto)
    ctodf = merge(ctodf, nodetissue)
    rownames(ctodf) = ctodf$gw.node
} else {
    diagcto = apply(lpmat.group, 1, which.max)
    ctodf = data.frame(node=names(diagcto), GROUP=odf$GROUP[diagcto], COLOR=odf$COLOR[diagcto])
    ctodf = unique(ctodf)
    rownames(ctodf) = ctodf$node
}


suid = "30643256 - Neuroticism"
suid = "29892015 - Atrial fibrillation"
suid = "30061737 - Atrial fibrillation"
suid = "30013184 - Allergic rhinitis"
suid = "30038396 - Educational attainment (MTAG)"
suid = "29212778 - Coronary artery disease"


klinks = filter(ddf, G1 == suid | G2 == suid)
knodes = sort(unique(c(as.character(klinks$G1), 
                       as.character(klinks$G2))))

lbs = nodes
lbs[!(lbs %in% knodes)] = ''
lbs = sapply(sub("^[0-9]* - ", "", lbs), width=20, split.text)

# Simple network: just the links/points:
sdf = ddf
net <- graph_from_data_frame(d=sdf, vertices=nodes, directed=F) 
vcol = as.character(ctodf[nodes,'COLOR'])
vcol[is.na(vcol)] = 'black'
ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
V(net)$size = 2
V(net)$label = lbs
V(net)$label.cex = .25
V(net)$label.color = rgb(0,0,0,.8)
V(net)$color = vcol
V(net)$frame.color <- 'black' # vcol
V(net)$frame.color <- NA
V(net)$pch = 19
E(net)$color = ecol 
if (by.zscore) {
    E(net)$width = sdf$zscore /5
E(net)$weight = sdf$zscore / 10
} else {
    E(net)$width = sdf$measure * 5
    E(net)$weight = sdf$measure / 10
}
set.seed(2)
l <- layout_with_fr(net, grid='nogrid') # Usually best
# l <- layout_components(net, layout=layout_with_fr, grid='nogrid') # decent alternative

png(paste0('~/', measure, '_network.png'),res=450,units='in',width=6,height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(net, layout=l)
dev.off()




