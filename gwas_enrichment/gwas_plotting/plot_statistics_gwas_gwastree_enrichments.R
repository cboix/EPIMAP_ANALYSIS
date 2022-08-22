#!/usr/bin/R
# -----------------------------------------
# Plot the gwas tree enrichment statistics:
# GWAS-centric analyses/plots:
# -----------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'great_enrichment/auxfunctions_GREAT.R'))
library(igraph)
library(flexclust)

# Arguments for loading data:
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions
use.adj = TRUE
# Cutoff at 1% with one FDR cutoff:
use.strict = FALSE
use.onecutoff = TRUE

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only, use.adj, use.strict, use.onecutoff) }
source(paste0(bindir, 'load_statistics_gwastree_enrichments.R'))

# Clean up:
rm(dflist)
rm(cdll)
rm(enhdf)
rm(regwide)
gc()

use.allgwas = TRUE
if (use.allgwas){
    # All gwas with at least 20k individuals:
    keep.cols = guid
    ksuff = 'allgwas'
} else {
    # Top sampsize traits with at least 20k individuals:
    keep.cols = kuid
    ksuff = 'traits'
}

use.cancersecondary = TRUE
if(use.cancersecondary){
    print("Marking cancer samples as their secondary tissue")
    imgpref = paste0(imgpref, 'cc2_')
    extpref = paste0(extpref, 'cc2_')

    nodetissue = unique(nodemeta[nodemeta$type == 'second', c('node','total','maxgroup')])
    names(nodetissue)[3] = 'GROUP'
    nodetissue =merge(nodetissue, rdcol, all.x=TRUE)
    nodetissue$COLOR[is.na(nodetissue$COLOR)] = 'grey80'  # Add color for Multiple
    nodetissue$category[is.na(nodetissue$category)] = 'Other'  # Add color for Multiple
    nodetissue.stats = aggregate(node ~ GROUP, nodetissue, length) # Only 191 assigned to 1

    ntmat = 1 * t(sapply(1:NN, function(x){odf$GROUP %in% nodetissue$GROUP[nodetissue$node == x]}))
    colnames(ntmat) = odf$GROUP
}

# Investigate diff:
serv.uids = scan('kept_uids_951_server.txt','c', sep="\n")
sum(keptgwas %in% serv.uids) # Fully in this
sum(serv.uids %in% keptgwas)
notmatched = serv.uids[!(serv.uids %in% keptgwas)]
# unmatched because of corrected sample size! 
# Need to transfer the new dataset to the server + re-sub all counts/ rer-run if correction will be affected
nsdf = gwssdf[gwssdf$uid %in% notmatched,]

# Export the gwas information:
ssdf = unique(gwssdf[gwssdf$uid %in% guid,c('uid','sampsize','pubMedID','trait','initSample','pubDate', 'pValue')])
sgdf = unique(gwdf[gwdf$uid %in% guid,c('chrom','chromStart','trait','pubMedID','pValue')])
save(ssdf, sgdf, file='kept_gwas_stats.Rda')


# -----------------------------------------------------
# Image of the kept GWAS, clustered by node occurrence:
# -----------------------------------------------------
redmat = all.regmat[keep.cols,]
redmat[redmat < 3] = 0
Zmat = 1 * (redmat >= 3)
leafmat = redmat %*% tmat
Zmat = 1 * (redmat >= 3) %*% ntmat

Zmat = sweep(Zmat, 1, apply(Zmat, 1, sum), '/')
sapply(c(0.5, 0.75, 1), function(x){sum(apply(Zmat, 1, max) >= x)})
# hist(apply(Zmat, 1, max), 20, col='grey')

rn = names(which(apply(Zmat, 1, max) < .50))
irn = names(which(apply(Zmat, 1, max) >= .50))
summary(keptgw$sampsize[keptgw$uid %in% rn])
summary(keptgw$sampsize[keptgw$uid %in% irn])
summary(apply(Zmat[rn,] > 0, 1, sum))

prn = names(which(apply(Zmat > 0, 1, sum) >= 34 * .25))
mean(apply(Zmat[prn,] > 0,1,sum))

mrn = rn[!(rn %in% prn)]
mean(apply(Zmat[mrn,] > 0,1,sum))

# Label mono-factorial, multi-factorial, poly-factorial.
# Give cutoffs:

# ONLY 1 tissue:
write.table(names(which(apply(Zmat, 1, max) == 1)), paste0(extpref, 'fdrp1_names_only_1.tsv'),
            quote=F, row.names=F, col.names=F, sep="\t")
# Above 75%:
write.table(names(which(apply(Zmat, 1, max) >= .75)), paste0(extpref, 'fdrp1_names_75over.tsv'),
            quote=F, row.names=F, col.names=F, sep="\t")
# Above 50%
write.table(names(which(apply(Zmat, 1, max) >= .50)),paste0(extpref, 'fdrp1_names_50over.tsv'),
            quote=F, row.names=F, col.names=F, sep="\t")
# Above 50% of GROUPS (poly-factorial)
write.table(names(which(apply(Zmat>0, 1, sum) >= 34 * .25)), paste0(extpref, 'fdrp1_names_9overgroup.tsv'),
            quote=F, row.names=F, col.names=F, sep="\t")
write.table(names(which(apply(Zmat>0, 1, sum) >= 34 * .5)), paste0(extpref, 'fdrp1_names_17overgroup.tsv'),
            quote=F, row.names=F, col.names=F, sep="\t")

sort(Zmat["29212778 - Coronary artery disease",] * 100, decreasing=T)

# Turn into tissue matrix - either normalize or reduce to hit:
# Ztmat = Zmat %*% ntmat
Ztmat = redmat %*% ntmat
Ztmat = sweep(Ztmat, 1, apply(Ztmat, 1, sum), '/')

# Statistics of single tissue/multi-tissue:

groupord = odf$GROUP
Ztmat = Ztmat[,groupord] # Reverse for plotting
# # Diagonalize:
dt = dist(Ztmat, 'ejaccard')
dt = dist(Ztmat, 'cosine')
ht <- hclust(dt, method='ward.D')
ht <- hclust(dt, method='complete')
cocl <- order.optimal(dt, ht$merge)$order
rn <- names(cocl)[cocl]
Ztmat = Ztmat[rn,]

ord = order(apply(Ztmat, 1, which.max), decreasing=TRUE)
Ztmat = Ztmat[ord,]
cto = apply(Ztmat, 1, which.max)
diagcto = cto
breaks = calc.breaks.acut(cto)
blabels = colnames(Ztmat)[unique(cto)]
bcol = odf$COLOR[unique(cto)]
ends = (c(0, nrow(Ztmat)) - .5) / (nrow(Ztmat) - 1)
bloc = (c(ends[1], breaks) + c(breaks, ends[2])) / 2
# Break by category diff:
colat = which(diff(as.numeric(odf$category)) == 1)
colbreaks = (colat - .5) / (ncol(Ztmat) - 1)

# Flip for plotting:
Ztmat = t(Ztmat)


# TODO: 1. Plot sections (3 col layout)
# TODO: 2. Output enrichments for people to use.

png(paste0(imgpref, 'gwas_red_clust_', ksuff, '.png'),res=450,units='in',width=5,height=18)
par(mar=c(.5,15,3,.5))
par(yaxs='i')
par(xaxs='i')
image(Ztmat, axes=F, col=col1)
yt = seq(0, 1, length.out=ncol(Ztmat))
xt = seq(0, 1, length.out=nrow(Ztmat))
yaxlab=colnames(Ztmat)
text(y=yt, x=par()$usr[1] - 0.01*(diff(par()$usr[1:2])),
     labels=colnames(Ztmat), srt=0, adj=1, xpd=TRUE,cex=1 / ncol(Ztmat) * 100)
text(x=xt, y=par()$usr[4] + 0.002*(diff(par()$usr[3:4])),
     labels=rownames(Ztmat), srt=90, adj=0, xpd=TRUE, cex=.35)
abline(h=breaks, lty='dotted')
abline(v=colbreaks, lty='dotted', lwd=.5)
box()
par(xpd=TRUE)
text(y=bloc, x=.02, adj=0,
     labels=blabels, col=bcol, cex=.75)
dev.off()

# Replace:
repl = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
              "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
colnames(Ztmat)[colnames(Ztmat) == repl] = "26974007 - Chronic inflammatory diseases (pleiotropy)"

plot_sub_Ztmat = function(pltid){
    # Make breaks, labels:
    szmat = Ztmat[,pltid]
    sbreaks = calc.breaks.acut(cto[pltid])
    sblabels = rownames(szmat)[unique(cto[pltid])]
    sbcol = odf$COLOR[unique(cto[pltid])]
    sends = (c(0, ncol(szmat)) - .5) / (ncol(szmat) - 1)
    sbloc = (c(sends[1], sbreaks) + c(sbreaks, sends[2])) / 2
    # Break by category diff:
    colat = which(diff(as.numeric(odf$category)) == 1)
    colbreaks = (colat - .5) / (nrow(szmat) - 1)
    # 
    par(mar=c(0.25,11,2.5,0.25))
    par(yaxs='i')
    par(xaxs='i')
    image(szmat, axes=F, col=col1, useRaster=TRUE)
    yt = seq(0, 1, length.out=ncol(szmat))
    xt = seq(0, 1, length.out=nrow(szmat))
    yaxlab=colnames(szmat)
    text(y=yt, x=par()$usr[1] - 0.01*(diff(par()$usr[1:2])),
         # labels=colnames(szmat), srt=0, adj=1, xpd=TRUE,cex=.5 / ncol(szmat) * 100)
         labels=colnames(szmat), srt=0, adj=1, xpd=TRUE,cex=.27)
    text(x=xt, y=par()$usr[4] + 0.002*(diff(par()$usr[3:4])),
         labels=rownames(szmat), srt=90, adj=0, xpd=TRUE, cex=.3)
    abline(h=sbreaks, lty='dashed', lwd=.25)
    abline(v=colbreaks, lty='dashed', lwd=.25)
    box(lwd=.5)
    par(xpd=TRUE)
    box.pad = 0.008
    xx = space.1d(sbloc, box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
    xx = space.1d(xx, box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
    text(y=xx, x=.02, adj=0, labels=sblabels, col=sbcol, cex=.35)
    par(xpd=FALSE)
}


pdf(paste0(imgpref, 'gwas_red_clust_', ksuff, '_1x3.pdf'),width=8,height=5)
layout(matrix(1:3, 1,3))
midSR = round(ncol(Ztmat) / 3)
plot_sub_Ztmat((2*midSR+1):ncol(Ztmat))
plot_sub_Ztmat((midSR+1):(2*midSR))
plot_sub_Ztmat(1:midSR)
dev.off()

pdf(paste0(imgpref, 'gwas_red_clust_', ksuff, '_1x2.pdf'),width=8,height=10.5)
layout(matrix(1:2, 1,2))
midSR = round(ncol(Ztmat) / 2)
plot_sub_Ztmat((midSR+1):ncol(Ztmat))
plot_sub_Ztmat(1:midSR)
dev.off()

Znam = colnames(Ztmat)
save(Znam, file=paste0(extpref, 'kept_', ksuff, '_ordered', cutsuff, suffix))

# ------------------------------
# Reduce, showing with wordlets:
# ------------------------------
NBIN = 2
wordlets <- get_summary_terms(colnames(Ztmat), NBIN)
wordlets_pos <- seq(1, ncol(Ztmat), length.out=length(wordlets)+1);
wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;

png(paste0(imgpref, 'gwas_red_clust_wdlt_', ksuff,'.png'),res=450,units='in',width=5,height=18)
par(mar=c(.5,15,3,.5))
par(yaxs='i')
par(xaxs='i')
image(Ztmat, axes=F, col=col1)
yt = seq(0, 1, length.out=ncol(Ztmat))
xt = seq(0, 1, length.out=nrow(Ztmat))
yaxlab=colnames(Ztmat)
text(x=xt, y=par()$usr[4] + 0.002*(diff(par()$usr[3:4])),
     labels=rownames(Ztmat), srt=90, adj=0, xpd=TRUE, cex=.35)
# text(y=yt, x=par()$usr[1] - 0.01*(diff(par()$usr[1:2])),
#      labels=colnames(Ztmat), srt=0, adj=1, xpd=TRUE,cex=1 / ncol(Ztmat) * 100)
text(y=wordlets_pos / ncol(Ztmat),
     x=par()$usr[1]-0.005*(par()$usr[2]-par()$usr[1]), 
     labels=wordlets, srt=0, adj=1, xpd=TRUE, 
     cex=1 / ncol(Ztmat) * 100 * 2)
abline(h=breaks, lty='dotted')
abline(v=colbreaks, lty='dotted', lwd=.5)
box()
par(xpd=TRUE)
text(y=bloc, x=.02, adj=0,
     labels=blabels, col=bcol, cex=.75)
dev.off()


# -------------------------
# Clustered non-singletons:
# -------------------------
# Also cluster traits into groups: (maybe groups like adipose + liver + muscle exist)
maxcut = 0.6
mind = which(apply(Ztmat, 2, max) <= maxcut)
Ctmat = Ztmat[, mind]
CZtmat = 1 * (Ctmat > 0.25)

dt = dist(t(CZtmat), 'Jaccard')
ht <- hclust(dt, method='ward.D')
# ht <- hclust(dt, method='complete')
cocl <- order.optimal(dt, ht$merge)$order
rn <- names(cocl)[cocl]
Ctmat = Ctmat[,rn]
CZtmat = CZtmat[,rn]
ord = order(apply(sweep(CZtmat, 1, 1:nrow(CZtmat), '*'), 2, mean), decreasing=T)
Ctmat = Ctmat[,ord]
CZtmat = CZtmat[,ord]
cto = apply(sweep(CZtmat, 1, 1:nrow(CZtmat), '*'), 2, mean)
breaks = calc.breaks.acut(round(cto, 1))

png(paste0(imgpref, 'gwas_red_clust_nonsingle_', ksuff, '.png'),res=450,units='in',width=5,height=18)
par(mar=c(.5,15,3,.5))
par(yaxs='i')
par(xaxs='i')
image(Ctmat, axes=F, col=col1)
yt = seq(0, 1, length.out=ncol(Ctmat))
xt = seq(0, 1, length.out=nrow(Ctmat))
yaxlab=colnames(Ctmat)
text(y=yt, x=par()$usr[1] - 0.01*(diff(par()$usr[1:2])),
     labels=colnames(Ctmat), srt=0, adj=1, xpd=TRUE,cex=1 / ncol(Ztmat) * 100)
text(x=xt, y=par()$usr[4] + 0.002*(diff(par()$usr[3:4])),
     labels=rownames(Ctmat), srt=90, adj=0, xpd=TRUE, cex=.35)
abline(h=breaks, lty='dotted')
abline(v=colbreaks, lty='dotted', lwd=.5)
box()
# par(xpd=TRUE)
# text(y=bloc, x=.02, adj=0, labels=blabels, col=bcol, cex=.75)
dev.off()


# ----------------------
# Clustered by KCCA GWAS matrix:
# ----------------------
fam <- 'ejaccard'
NCENTERS = 50
cl1 <- kcca(t(Ztmat), k=NCENTERS, family=kccaFamily(fam))
summary(cl1)

# Get kcca cluster attributes:
acl <- attributes(cl1)
dcl <- dist(acl$centers,'ejaccard')
hcl <- hclust(dcl)
# Order by max:
acl$centers
cocl = order(apply(acl$centers, 1, which.max), decreasing=TRUE)
cls <- ordered(acl$cluster,levels=cocl)
ord <- order(cls) # will order according to factor
cto = as.numeric(cls[ord])
breaks = calc.breaks.acut(cto)

# Order by diag
Ztmat = Ztmat[,ord]

png(paste0(imgpref, 'gwas_red_kcca_clust_', ksuff, '.png'),res=450,units='in',width=9,height=18)
par(mar=c(.5,15,7,.5))
par(yaxs='i')
par(xaxs='i')
image(Ztmat, axes=F, col=col1)
yt = seq(0, 1, length.out=ncol(Ztmat))
xt = seq(0, 1, length.out=nrow(Ztmat))
yaxlab=colnames(Ztmat)
text(y=yt, x=par()$usr[1] - 0.01*(diff(par()$usr[1:2])),
     labels=colnames(Ztmat), srt=0, adj=1, xpd=TRUE,cex=1 / ncol(Ztmat) * 100)
text(x=xt, y=par()$usr[4] + 0.005*(diff(par()$usr[3:4])),
     labels=rownames(Ztmat), srt=90, adj=0, xpd=TRUE, cex=1)
abline(h=breaks, lty='dotted')
abline(v=colbreaks, lty='dotted', lwd=.5)
box()
par(xpd=TRUE)
text(y=bloc, x=.01, adj=0,
     labels=blabels, col=bcol, cex=.75)
dev.off()





# -------------------------------------
# GWAS network: 
# Color is maximal of mult
# Adjacency is cosine distance
# -------------------------------------
# Use the normalized node matrix:
# TODO: do the sharing/colors directly on the node x GWAS matrix?
# binary.mat=TRUE # Old figs not binary.
binary.mat=FALSE # Old figs not binary.
all.nodes=FALSE
# metric = 'jaccard'
metric = 'cosine'
npref = paste('gwas', metric, 'network', sep="_")
if (binary.mat){
    Ztmat = 1 * (all.regmat[guid,] > 0) %*% ntmat
    npref = paste0(npref, '_binary')
} else {
    Ztmat = all.regmat[guid,] %*% ntmat
}
if (metric == 'jaccard'){ 
    Ztmat[Ztmat > 0.05] = 1 
    cutoff = 0.49
} else {
    cutoff = 0.75
    # cutoff = 0.7
    # cutoff = 0.65
}
npref = paste0(npref, '_', sub("\\.", "_", cutoff))

dt = dist(Ztmat, metric)
dmat = 1 - as.matrix(dt)
rn = rownames(dmat)
dmat = data.frame(dmat)
colnames(dmat) = rn
rownames(dmat) = rn
dmat$T1 = rownames(dmat)
ddf  = gather(dmat, T2, sim, -T1)
ddf = ddf[ddf$sim >= cutoff,]
ddf = ddf[ddf$T1 != ddf$T2,]
if (all.nodes){
    nodes = sort(guid)
    npref = paste0(npref, '_allnodes')
} else {
    nodes = sort(unique(c(ddf$T1, ddf$T2)))
}
# Remove edges in opposite direction
ddf$T1 = factor(ddf$T1, levels=nodes)
ddf$T2 = factor(ddf$T2, levels=nodes)
ddf = ddf[as.numeric(ddf$T1) < as.numeric(ddf$T2),]
# Kept pct:
dim(ddf)[1]
dim(ddf)[1] / (dim(dmat)[1]^2)
m = apply(ddf, 1, function(x){which.max(Ztmat[x[1],] * Ztmat[x[2],])})
ddf$COLOR = odf$COLOR[m]
diagcto = apply(Ztmat, 1, which.max)
ctodf = data.frame(node=names(diagcto), GROUP=odf$GROUP[diagcto], COLOR=odf$COLOR[diagcto])
rownames(ctodf) = ctodf$node

sum(ctodf[as.character(ddf$T1), 'COLOR'] ==
    ctodf[as.character(ddf$T2), 'COLOR']) /nrow(ddf)

# Simple network: just the links/points:
sdf = ddf
net <- graph_from_data_frame(d=sdf, vertices=nodes, directed=F) 
vcol = as.character(ctodf[nodes,'COLOR'])
vcol[is.na(vcol)] = 'black'
ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
V(net)$size = 2
V(net)$label = ""
V(net)$color = vcol
V(net)$frame.color <- 'black' # vcol
V(net)$frame.color <- NA
V(net)$pch = 19
E(net)$color = ecol 
elty = rep('dotted', length(sdf$sim))
elty[sdf$sim >= .85] = 'dashed'
elty[sdf$sim >= .95] = 'solid'
E(net)$lty = elty
ewidth = (sdf$sim >= .95) * .4 + (sdf$sim >= .85) * 0.4 + (sdf$sim >= 0.75) * 0.4 + 0.4
E(net)$width = sdf$sim
# E(net)$weight = sdf$sim * .5
E(net)$weight = sdf$sim * .075
# E(net)$weight = sdf$sim * .5
# set.seed(2)
set.seed(8)
l <- layout_with_fr(net, grid='nogrid') # Usually best

# png(paste0(imgpref, npref, '.png'),res=450,units='in',width=6,height=6)
pdf(paste0(imgpref, npref, '.pdf'), width=6, height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(net, layout=l, curved=F)
# plot(net, layout=l, edge.curved=seq(-0.5, 0.5, length = ecount(net)))
dev.off()

# Pie chart experimental
pie.values = lapply(nodes, function(x){ Ztmat[x,] })
maxp = apply(all.regmat[nodes,], 1, max)
bkcut = c(3, 10, 25, 50, 1000)
fmaxp = cut(maxp, breaks=bkcut, right=FALSE)
maxv = (as.numeric(fmaxp) + 1) / 2

pdf(paste0(imgpref, npref, '_pie_linetypes_discretesizepts.pdf'),width=6,height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(net, vertex.shape="pie", vertex.pie=pie.values, layout=l,
     vertex.pie.color=list(odf$COLOR),
     vertex.pie.border=vcol,
     vertex.pie.lty=1,
     vertex.size=maxv)
legend('bottomright', 
       legend=c('max -log10p >= 3', 'max -log10p >= 10', 'max -log10p >= 25', 'max -log10p >= 50', 
                'Cosine Sim. >= 0.75', 'Cosine Sim. >= 0.85', 'Cosine Sim. >= 0.95'),
       col='black', lty=c(NA, NA, NA, NA, 'dotted','dashed','solid'),
       cex=.4, pch=c(19, 19,19,19, NA, NA, NA), pt.cex=c(1, 1.5, 2, 2.5, .75, .75, .75) / 4, 
       inset=0, box.col=NA)
legend('bottomleft', legend=odf$GROUP, col=odf$COLOR, ncol=5,
       cex=.4, pch=15, pt.cex=1, inset=0, box.col=NA)
dev.off()


maxv = (log10(maxp) + 0.5) * 1.25
csizes = sapply(bkcut[1:4], function(x){ min(maxv[which(maxp >= x)])})

pdf(paste0(imgpref, npref, '_pie_linetypes_continsizepts.pdf'),width=6,height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(net, vertex.shape="pie", vertex.pie=pie.values, layout=l,
     vertex.pie.color=list(odf$COLOR),
     vertex.pie.border=vcol,
     vertex.pie.lty=1,
     vertex.size=maxv)
legend('bottomright', 
       legend=c('max -log10p >= 3', 'max -log10p >= 10', 'max -log10p >= 25', 'max -log10p >= 50', 
                'Cosine Sim. >= 0.75', 'Cosine Sim. >= 0.85', 'Cosine Sim. >= 0.95'),
       col='black', lty=c(NA, NA, NA, NA, 'dotted','dashed','solid'),
       cex=.4, pch=c(19, 19,19,19, NA, NA, NA), pt.cex=c(csizes, NA, NA, NA)/ 4, 
       inset=0, box.col=NA)
legend('bottomleft', legend=odf$GROUP, col=odf$COLOR, ncol=5,
       cex=.4, pch=19, pt.cex=.4, inset=0, box.col=NA)
dev.off()


wcut = c(.75, .85, .95, 1.05)
fwidth = cut(sdf$sim, breaks=wcut, right=FALSE)
ewidth = as.numeric(fwidth) * .5
E(net)$lty = 'solid'
E(net)$width = ewidth
esizes = sort(unique(ewidth))

pdf(paste0(imgpref, npref, '_pie_linewidths_continsizepts.pdf'),width=6,height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(net, vertex.shape="pie", vertex.pie=pie.values, layout=l,
     vertex.pie.color=list(odf$COLOR),
     vertex.pie.border=vcol,
     vertex.pie.lty=1,
     vertex.size=maxv)
legend('bottomright', 
       legend=c('max -log10p >= 3', 'max -log10p >= 10', 'max -log10p >= 25', 'max -log10p >= 50', 
                'Cosine Sim. >= 0.75', 'Cosine Sim. >= 0.85', 'Cosine Sim. >= 0.95'),
       col='black', lty=c(NA, NA, NA, NA, 'solid','solid','solid'), lwd=c(rep(NA, 4), esizes),
       cex=.4, pch=c(19, 19,19,19, NA, NA, NA), pt.cex=c(csizes, NA, NA, NA)/ 4, 
       inset=0, box.col=NA)
legend('bottomleft', legend=odf$GROUP, col=odf$COLOR, ncol=5,
       cex=.4, pch=19, pt.cex=.4, inset=0, box.col=NA)
dev.off()



# -------------------------
# Network with some labels:
# -------------------------
# TODO: keep labels on the ones in the middle of things
# Labeling: 
# 1. Keep highest pvalue nodes per trait where the trait is top trait:
# Nodes w/ trait top:
Ztmat = all.regmat[guid,] %*% ntmat
nZtmat = sweep(Ztmat, 1, apply(Ztmat,1,sum), '/')
tp.df = data.frame(uid = rownames(Ztmat), 
                   pagg = apply(Ztmat, 1, max),
                   GROUP=colnames(Ztmat)[apply(Ztmat, 1, which.max)], 
                   frac = apply(nZtmat, 1, max))
tp.df = tp.df[tp.df$frac > 0.50,]
tp.df = tp.df[tp.df$pagg > 10,]
take = 2
tuid = c()
for (g in odf$GROUP){
    subdf = tp.df[tp.df$GROUP == g,]
    subdf = subdf[order(subdf$pagg, decreasing=T), ]
    tuid = c(tuid, head(as.character(subdf$uid), take))
}
tp.df = tp.df[tp.df$uid %in% tuid,]
tp.df$trait = sapply(sub("^[0-9]* - ", "",tp.df$uid), width=40, split.text)
hubnodes = as.character(tp.df$uid)

# Find mixed nodes:
wind = which(apply(Ztmat,1, max) > 10)
wsec =apply(nZtmat[wind,], 1, function(x){sort(x, decreasing=T)[2]})
linkhits = names(which(wsec > .33))
# linkhits = linkhits[nchar(linkhits)]

# TODO: want to space better + color the relevant GWAS

# 3. Label near center-mass for each trait
# Get the centermass:
lablocs = c()
glabs = c()
for (sgroup in odf$GROUP){
    ind = which(ctodf[nodes, 'GROUP'] == sgroup)
    if (length(ind) > 1){
        lablocs = rbind(lablocs, apply(l[ind,], 2, median))
        glabs = c(glabs, sgroup)
    }
    # NOTE: Other is poor qual. want to map to representative.
}
# Scale locations:
lrange = apply(l, 2, range)
lablocs = sweep(lablocs, 2, lrange[1,], '-')
lablocs = sweep(lablocs, 2, lrange[2,] - lrange[1,], '/') * 2 - 1
# Move labels outwards:
lablocs = lablocs * 1.1
lablocs[lablocs > 1] = 1
lablocs[lablocs < -1] = -1

# Color hubs and mid differently.
# kind = c(which(nodes %in% hubnodes), which(nodes %in% lindf$uid))
# kind = c(which(nodes %in% hubnodes))
kind = c(which(nodes %in% hubnodes), which(nodes %in% linkhits))
lnet = net
traits = sub("^[0-9][0-9]* - ", "", nodes)
labs = rep("", length(nodes))
labs[kind] = traits[kind]
labs = sapply(labs, width=50,split.text)
# V(lnet)$label = labs
V(lnet)$label = NA
V(lnet)$label.cex = .4
V(lnet)$label.color = rgb(0,0,0,.5)

# Use space 1x in both directions using text width and height.
# text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
#      srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))


png(paste0(imgpref, npref, '_labels.png'),res=450,units='in',width=6,height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(lnet, layout=l, curved=F)
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
dev.off()

png(paste0(imgpref, npref, '_nolabels_mainonly.png'), res=450, units='in', width=6,height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(lnet, layout=l, curved=F)
dev.off()

pdf(paste0(imgpref, npref, '_onlylabels_mainonly.pdf'), width=6,height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(1,1, xlim=c(-1,1), ylim=c(-1,1), type='n')
# plot(lnet, layout=l, curved=F)
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
dev.off()

png(paste0(imgpref, npref, '_pie_labels.png'),res=450,units='in',width=6,height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(lnet, vertex.shape="pie", vertex.pie=pie.values, layout=l,
     vertex.pie.color=list(odf$COLOR),
     vertex.pie.border=vcol,
     vertex.pie.lty=1,
     vertex.size=maxv)
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
dev.off()



traits = sub("^[0-9][0-9]* - ", "", nodes)
labs = traits
labs = sapply(labs, width=50,split.text)
V(lnet)$label = labs
V(lnet)$label.cex = .25

# png(paste0(imgpref, npref, '_pie_labels_traits.png'),res=450,units='in',width=6,height=6)
srt.maxp = sort(maxp, decreasing=T)
t50 = head(srt.maxp,50)[50]
t100 = head(srt.maxp,100)[100]
t200 = head(srt.maxp,200)[200]
t300 = head(srt.maxp,300)[300]
vlabcol = ifelse(maxp > t50, 'black', ifelse(maxp > t100, 'grey15', ifelse(maxp > t200, 'grey30', ifelse(maxp > t300, 'grey50','grey70'))))
png(paste0(imgpref, npref, '_pie_labels_traits.png'),res=450,units='in',width=10,height=10)
sp = 0.1
par(mar = rep(sp,4))
plot(lnet, vertex.shape="pie", vertex.pie=pie.values, layout=l,
     vertex.pie.color=list(odf$COLOR),
     vertex.pie.border=vcol,
     vertex.pie.lty=1,
     vertex.label.color=vlabcol,
     vertex.size=maxv)
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
dev.off()

# png(paste0(imgpref, npref, '_pie_labels_traits.png'),res=450,units='in',width=6,height=6)
pdf(paste0(imgpref, npref, '_pie_labels_traits.pdf'), width=10, height=10)
sp = 0.1
par(mar = rep(sp,4))
plot(lnet, vertex.shape="pie", vertex.pie=pie.values, layout=l,
     vertex.pie.color=list(odf$COLOR),
     vertex.pie.border=vcol,
     vertex.pie.lty=1,
     vertex.label.color=vlabcol,
     vertex.size=maxv)
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
dev.off()


# REPEL LABELS:
source(paste0(bindir, 'auxiliary_function_general_repel.R'))
TOPBETWEEN = 100
btlnet = betweenness(lnet)
topbt = names(sort(btlnet, decreasing=T)[1:TOPBETWEEN])
# TOP PVALS
# Add
TOPPVALS = 250
Zmat = sweep(Ztmat, 1, apply(Ztmat, 1, sum), '/')
Zmat = Zmat[nodes,]
zind = rownames(Zmat)[!(rownames(Zmat) %in% topbt)]
subZmat = Zmat[zind,]

# TOP p-vals:
toppval = names(head(sort(apply(subZmat, 1, max)), TOPPVALS))
# TOP BETWEENNESS
traits = sub("^[0-9][0-9]* - ", "", nodes)
# Pull only top 50 by betweenness centrality:
topind = which(nodes %in% c(topbt, toppval))
traits[-topind] = ""
print(length(topind))
top.maxp = maxp[topind]
top.vlabcol = vlabcol[topind]

# Network params:
V(lnet)$label = NA
l2 = l
l2 = sweep(l2, 2, lrange[1,], '-')
l2 = sweep(l2, 2, lrange[2,] - lrange[1,], '/') * 2 - 1

numpref = paste0('btwn', TOPBETWEEN, '_pval', TOPPVALS)
# png(paste0(imgpref, npref, '_pie_labels_traits.png'),res=450,units='in',width=6,height=6)
pdf(paste0(imgpref, npref, '_pie_labels_traits_repel_', numpref, '.pdf'), width=10, height=10)
sp = 2
par(mar = rep(sp,4))
plot(lnet, vertex.shape="pie", vertex.pie=pie.values, layout=l,
     vertex.pie.color=list(odf$COLOR), vertex.pie.border=vcol,
     vertex.pie.lty=1, vertex.size=maxv)
# Place labels:
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
print(par()$usr)
# Repel points:
lbcex=0.25
rdf = general_repel_text(x=l2[,1], y=l2[,2], 
                         xlim=par()$usr[1:2] * 1.25, ylim=par()$usr[3:4] * 1.25,
                         hjust=.5, vjust=.5, seed=1, max.iter=5000,
                         labels=traits, cex=lbcex, pt.cex=.25)
text(x=rdf$x, y=rdf$y, labels=rdf$lab,
     srt=0, adj=0, xpd=TRUE, cex=lbcex, col=top.vlabcol)
segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25, col=top.vlabcol)
dev.off()


# png(paste0(imgpref, npref, '_pie_labels_traits.png'),res=450,units='in',width=6,height=6)
pdf(paste0(imgpref, npref, '_labels_traits_repel_', numpref, '.pdf'), width=10, height=10)
sp = 2
par(mar = rep(sp,4))
plot(lnet, layout=l, curved=F)
# Place labels:
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
# Repel points:
lbcex=0.25
rdf = general_repel_text(x=l2[,1], y=l2[,2], hjust=.5, vjust=.5, seed=1,
                         labels=traits, cex=lbcex, pt.cex=.25)
text(x=rdf$x, y=rdf$y, labels=rdf$lab,
     srt=0, adj=0, xpd=TRUE, cex=lbcex, col='black')
segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25, col='grey')
dev.off()


# ALL labels:
# png(paste0(imgpref, npref, '_pie_labels_traits.png'),res=450,units='in',width=6,height=6)
pdf(paste0(imgpref, npref, '_pie_labels_traits_repel_allnonsingle.pdf'), width=10, height=10)
sp = 2
par(mar = rep(sp,4))
plot(lnet, vertex.shape="pie", vertex.pie=pie.values, layout=l,
     vertex.pie.color=list(odf$COLOR), vertex.pie.border=vcol,
     vertex.pie.lty=1, vertex.size=maxv)
# Place labels:
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
print(par()$usr)
# Repel points:
lbcex=0.25
rdf = general_repel_text(x=l2[,1], y=l2[,2], 
                         xlim=par()$usr[1:2] * 1.25, ylim=par()$usr[3:4] * 1.25,
                         hjust=.5, vjust=.5, seed=1, max.iter=5000,
                         labels=sub("^[0-9][0-9]* - ", "", nodes), cex=lbcex, pt.cex=.25)
text(x=rdf$x, y=rdf$y, labels=rdf$lab,
     srt=0, adj=0, xpd=TRUE, cex=lbcex, col=vlabcol)
segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25, col=vlabcol)
dev.off()


# ALL labels:
# png(paste0(imgpref, npref, '_pie_labels_traits.png'),res=450,units='in',width=6,height=6)
pdf(paste0(imgpref, npref, '_pie_labels_traits_repel_', numpref,'_remainder.pdf'), width=10, height=10)
sp = 2
par(mar = rep(sp,4))
plot(lnet, vertex.shape="pie", vertex.pie=pie.values, layout=l,
     vertex.pie.color=list(odf$COLOR), vertex.pie.border=vcol,
     vertex.pie.lty=1, vertex.size=maxv)
# Place labels:
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
print(par()$usr)
# Repel points:
rtraits = sub("^[0-9][0-9]* - ", "", nodes)
rtraits[topind] = ''
r.vlabcol = vlabcol[-topind]
lbcex=0.25
rdf = general_repel_text(x=l2[,1], y=l2[,2], 
                         xlim=par()$usr[1:2] * 1.25, ylim=par()$usr[3:4] * 1.25,
                         hjust=.5, vjust=.5, seed=1, max.iter=5000,
                         labels=rtraits, cex=lbcex, pt.cex=.25)
text(x=rdf$x, y=rdf$y, labels=rdf$lab,
     srt=0, adj=0, xpd=TRUE, cex=lbcex, col=r.vlabcol)
segments(rdf$x, rdf$y, rdf$x.orig, rdf$y.orig, lwd=.25, col=r.vlabcol)
dev.off()





# ----------------------------------
# Compare to the genetic similarity:
# ----------------------------------
# Rough 1st pass: bin SNPs into 10kb bins:
# Second pass should be - look at larger set
binsize = 10000
gwdf$bin = round(gwdf$chromStart / binsize, 0)
gwdf$sid = paste0(gwdf$chrom, "_", gwdf$bin)
# Unique hits:
dim(gwdf)
dim(unique(gwdf[,c('chrom','bin')])) # about halved
subgwdf = gwdf[gwdf$uid %in% keptgw$uid, ] 
dim(subgwdf)
dim(unique(subgwdf[,c('chrom','bin')])) # about halved

# Keep double snp:
subgwdf = gwdf[gwdf$uid %in% nodes,]
ksnp = aggregate(uid ~ sid, subgwdf, length)
ksnp = ksnp[ksnp$uid > 1,] 
dim(ksnp) # About 10k shared in the nodes set

# Genetics matrix:
gmat = matrix(0, nrow=length(nodes), ncol=length(nodes), dimnames=list(nodes, nodes))
for (i in 1:nrow(ksnp)){
    ssid = ksnp$sid[i]
    suids = as.character(subgwdf$uid[subgwdf$sid == ssid])
    gmat[t(combn(suids,2))] = gmat[t(combn(suids,2))] + 1
}

# Normalize the gmat by total pvalue:
nsnpdf = aggregate(pValue ~ uid, gwdf, length)
rownames(nsnpdf) = nsnpdf$uid
rs = nsnpdf[nodes, 'pValue']
rsm = matrix(rep(rs, length(nodes)), nrow=length(nodes),
             ncol=length(nodes), byrow=T, dimnames=list(nodes, nodes))
rsm = rsm + t(rsm)
nunion = rsm - gmat
gsim = gmat / nunion 
gsim = data.frame(gsim)

rn = rownames(gsim)
colnames(gsim) = rn
rownames(gsim) = rn
gsim$T1 = rownames(gsim)
gdf  = gather(gsim, T2, gensim, -T1)
cutoff = 0.05
gdf = gdf[gdf$gensim > cutoff,]
gdf = gdf[gdf$T1 != gdf$T2,]

# Remove edges in opposite direction
gdf$T1 = factor(gdf$T1, levels=nodes)
gdf$T2 = factor(gdf$T2, levels=nodes)
gdf = gdf[as.numeric(gdf$T1) < as.numeric(gdf$T2),]
# Kept pct:
dim(gdf)[1]
dim(gdf)[1] / (dim(gmat)[1]^2)
m = apply(gdf, 1, function(x){which.max(Ztmat[x[1],] * Ztmat[x[2],])})
gdf$COLOR = odf$COLOR[m]
ctodf = data.frame(node=names(diagcto), GROUP=odf$GROUP[diagcto], COLOR=odf$COLOR[diagcto])
rownames(ctodf) = ctodf$node

# Network from genetics only
sdf = gdf
gonet <- graph_from_data_frame(d=sdf, vertices=nodes, directed=F) 
ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
V(gonet)$size = 2
V(gonet)$label = ""
V(gonet)$frame.color = NA
V(gonet)$color = as.character(ctodf[nodes,'COLOR'])
V(gonet)$bg = as.character(ctodf[nodes,'COLOR'])
E(gonet)$color = ecol
E(gonet)$width = 1
l2 <- layout_with_fr(gonet)

# Network from genetics only + intersect with epigenetics
sdf = merge(ddf, gdf)
gnet <- graph_from_data_frame(d=sdf, vertices=nodes, directed=F) 
ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
V(gnet)$size = 2
V(gnet)$label = ""
V(gnet)$frame.color = NA
V(gnet)$color = as.character(ctodf[nodes,'COLOR'])
V(gnet)$bg = as.character(ctodf[nodes,'COLOR'])
E(gnet)$color = ecol
E(gnet)$width = 1

# Network from genetics only (only linked nodes)
sdf = gdf
genet.nodes = sort(unique(c(as.character(sdf$T1), 
                            as.character(sdf$T2))))
red.gonet <- graph_from_data_frame(d=sdf, vertices=genet.nodes, directed=F) 
ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
V(red.gonet)$size = 2
V(red.gonet)$label = ""
V(red.gonet)$frame.color = NA
V(red.gonet)$color = as.character(ctodf[genet.nodes,'COLOR'])
V(red.gonet)$bg = as.character(ctodf[genet.nodes,'COLOR'])
E(red.gonet)$color = ecol
E(red.gonet)$width = 1
l2 <- layout_with_fr(red.gonet)

# Simple network: just the links/points:
sdf = ddf[ddf$T1 %in% genet.nodes & ddf$T2 %in% genet.nodes,]
red.enet <- graph_from_data_frame(d=sdf, vertices=genet.nodes, directed=F) 
vcol = as.character(ctodf[genet.nodes,'COLOR'])
vcol[is.na(vcol)] = 'black'
ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
V(red.enet)$size = 2
V(red.enet)$label = ""
V(red.enet)$color = vcol
V(red.enet)$frame.color <- 'black' # vcol
V(red.enet)$frame.color <- NA
V(red.enet)$pch = 19
E(red.enet)$color = ecol 
E(red.enet)$width = sdf$sim
E(red.enet)$weight = sdf$sim * .5

png(paste0(imgpref, npref, '_genetics_comparison.png'),res=250,units='in',width=20,height=4)
layout(matrix(c(1:5), 1, 5), TRUE)
sp = 0.15
mars = c(sp,sp, 2, sp)
par(mar=mars)
plot(net, layout=l, curved=F)
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
mtext("Epigenetic overlap", side=3, line=0)
par(mar=mars)
plot(gnet, layout=l, curved=F)
mtext("Genetic epigenetic overlap", side=3, line=0)
par(mar=mars)
plot(gonet, layout=l, curved=F)
mtext("Any genetic overlap", side=3, line=0)
par(mar=mars)
plot(red.gonet, layout=l2, curved=F)
mtext("Any genetic overlap (layout genetics)", side=3, line=0)
par(mar=mars)
plot(red.enet, layout=l2, curved=F)
mtext("Epigenetic overlap (layout genetics)", side=3, line=0)
dev.off()

title.cex=0.5
png(paste0(imgpref, npref, '_genetics_comparison_2x2.png'),res=250,units='in',width=8,height=8)
layout(matrix(c(1:4), 2, 2), TRUE)
sp = 0.15
mars = c(sp,sp, 2, sp)
# par(mar=mars)
# plot(net, layout=l, curved=F)
# mtext("Epigenetic overlap", side=3, line=0)
par(mar=mars)
plot(gnet, layout=l, curved=F)
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
mtext("Epigenetic and Genetic Similarity", side=3, line=0.5, cex=title.cex)
mtext("Network Layout: Epigenetic Similarity", side=3, line=-.5, cex=0.8 * title.cex)
par(mar=mars)
plot(gonet, layout=l, curved=F)
mtext("Genetic Similarity", side=3, line=0.5, cex=title.cex)
mtext("Network Layout: Epigenetic Similarity", side=3, line=-.5, cex=0.8 * title.cex)
par(mar=mars)
plot(red.gonet, layout=l2, curved=F)
mtext("Genetic Similarity", side=3, line=0.5, cex=title.cex)
mtext("Network Layout: Genetic Similarity", side=3, line=-.5, cex=0.8 * title.cex)
par(mar=mars)
plot(red.enet, layout=l2, curved=F)
mtext("Epigenetic Similarity", side=3, line=0.5, cex=title.cex)
mtext("Network Layout: Genetic Similarity", side=3, line=-.5, cex=0.8 * title.cex)
dev.off()

png(paste0(imgpref, npref, '_genetics_comparison_panel.png'),res=450,units='in',width=8,height=4)
layout(matrix(c(1,1, 2:5), 2, 3), widths=c(4,2,2), heights=c(2,2), TRUE)
sp = 0.05
mars = c(sp,sp, 2, sp)
par(mar=mars)
plot(lnet, layout=l, curved=F)
mtext("Epigenetic overlap", side=3, line=0, cex=title.cex * 1.5)
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
par(mar=mars)
plot(gnet, layout=l, curved=F)
mtext("Genetic epigenetic overlap", side=3, line=0, cex=title.cex)
par(mar=mars)
plot(gonet, layout=l, curved=F)
mtext("Any genetic overlap", side=3, line=0, cex=title.cex)
par(mar=mars)
plot(red.gonet, layout=l2, curved=F)
mtext("Any genetic overlap (layout genetics)", side=3, line=0, cex=title.cex)
par(mar=mars)
plot(red.enet, layout=l2, curved=F)
mtext("Epigenetic overlap (layout genetics)", side=3, line=0, cex=title.cex)
dev.off()


png(paste0(imgpref, npref, '_genetics_comparison_panel_ptonly.png'),res=450,units='in',width=8,height=4)
layout(matrix(c(1,1, 2:5), 2, 3), widths=c(4,2,2), heights=c(2,2), TRUE)
sp = 0.05
title.cex=0.5
mars = c(sp,sp, 2, sp)
par(mar=mars)
plot(lnet, layout=l, curved=F)
# mtext("Epigenetic overlap", side=3, line=0, cex=title.cex * 1.5)
# text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
#      srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
par(mar=mars)
plot(gnet, layout=l, curved=F)
# mtext("Genetic epigenetic overlap", side=3, line=0, cex=title.cex)
par(mar=mars)
plot(gonet, layout=l, curved=F)
# mtext("Any genetic overlap", side=3, line=0, cex=title.cex)
par(mar=mars)
plot(red.gonet, layout=l2, curved=F)
# mtext("Any genetic overlap (layout genetics)", side=3, line=0, cex=title.cex)
par(mar=mars)
plot(red.enet, layout=l2, curved=F)
# mtext("Epigenetic overlap (layout genetics)", side=3, line=0, cex=title.cex)
dev.off()


pdf(paste0(imgpref, npref, '_genetics_comparison_panel_axonly.pdf'),width=8,height=4)
layout(matrix(c(1,1, 2:5), 2, 3), widths=c(4,2,2), heights=c(2,2), TRUE)
sp = 0.05
title.cex=0.75
mars = c(sp,sp, 2, sp)
par(mar=mars)
plot(1,1, xlim=c(-1,1), ylim=c(-1,1), type='n', axes=F)
mtext("Trait-Trait Similarity by Epigenetic Enrichment", side=3, line=0, cex=title.cex * 1.25)
text(x=lablocs[,1], y=lablocs[,2], labels=glabs, 
     srt=0, adj=0, xpd=TRUE,cex=.75, col=rgb(0,0,0,.5))
par(mar=mars)
plot(1,1, xlim=c(-1,1), ylim=c(-1,1), type='n', axes=F)
mtext("Epigenetic and Genetic Similarity", side=3, line=0.5, cex=title.cex)
mtext("Network Layout: Epigenetic Similarity", side=3, line=-.5, cex=0.8 * title.cex)
par(mar=mars)
plot(1,1, xlim=c(-1,1), ylim=c(-1,1), type='n', axes=F)
# mtext("Genetic Similarity\nLayout by Epigenetic Similarity", side=3, line=0, cex=title.cex)
mtext("Genetic Similarity", side=3, line=0.5, cex=title.cex)
mtext("Network Layout: Epigenetic Similarity", side=3, line=-.5, cex=0.8 * title.cex)
par(mar=mars)
plot(1,1, xlim=c(-1,1), ylim=c(-1,1), type='n', axes=F)
mtext("Genetic Similarity", side=3, line=0.5, cex=title.cex)
mtext("Network Layout: Genetic Similarity", side=3, line=-.5, cex=0.8 * title.cex)
par(mar=mars)
plot(1,1, xlim=c(-1,1), ylim=c(-1,1), type='n', axes=F)
mtext("Epigenetic Similarity", side=3, line=0.5, cex=title.cex)
mtext("Network Layout: Genetic Similarity", side=3, line=-.5, cex=0.8 * title.cex)
dev.off()






# TODO: Now list ones with genetic overlap but not epigenetic, and vice versa: 
# - basically going at question of pleiotropy / regulation.

# --------------------------
# Show these using heatmaps:
# --------------------------
colramp = viridis(100)
epimat = as.matrix(dmat[,-ncol(dmat)])
genmat = as.matrix(gsim[,-ncol(gsim)])

rn = rownames(epimat)
# Map gen to epi loc:
dt = as.dist(1 - epimat)
ht <- hclust(dt, method='ward.D')
cocl <- order.optimal(dt, ht$merge)$order
rn <- names(cocl)[cocl]

epimat = epimat[rn,rn]
grn = rn[rn %in% rownames(genmat)]
tmp = epimat
tmp[] = NA
tmp[grn, grn] = genmat[grn, grn]

diag(tmp) = NA
diag(epimat) = NA
vcol.f = factor(ctodf[rn,'COLOR'])
vcol = as.character(ctodf[rn,'COLOR'])

# Generate breaks:
NCLUST = 20
breaks = calc.breaks(ht, NCLUST, cocl)
acut <- cutree(ht, NCLUST)[cocl]

colramp = col1
# png(paste0(imgpref, 'gwas_genet_cosine_heatmap.png'),res=450,units='in',width=18,height=8)
pdf(paste0(imgpref, 'gwas_genet_cosine_heatmap.pdf'),width=12,height=6)
layout(matrix(c(1:4), 1, 4),  widths=c(.4,3,.1, 3), heights=3, TRUE)
rsp=5
sp=0.1
par(yaxs='i')
par(xaxs='i')
minor = 1
par(mar=c(sp,rsp,1,0))
image(t(as.matrix(as.numeric(vcol.f))), col=levels(vcol.f), useRaster=T, axes=F)
text(x=par()$usr[1] - 0.01 * diff(par()$usr[1:2]),
     y=seq(0,1,length.out=nrow(epimat)), 
     rownames(epimat), cex=0.1, adj=1, xpd=TRUE)
par(mar=c(sp,sp,1,sp))
# Color sides by the vcol
image(epimat, axes=F, col=c('white', col1), useRaster=TRUE)
mtext("Epigenetic enrichment similarity", side=3, cex=.75)
abline(v=breaks,lty=2,lwd=.5)
abline(h=breaks,lty=2,lwd=.5)
box(lwd=.5)
par(mar=c(sp,0,1,0))
image(t(as.matrix(as.numeric(vcol.f))), col=levels(vcol.f), useRaster=T, axes=F)
par(mar=c(sp,sp,1,sp))
image(1 * (tmp > 0), axes=F, col=c('white', col1), useRaster=TRUE)
abline(v=breaks,lty=2,lwd=.5)
abline(h=breaks,lty=2,lwd=.5)
box(lwd=.5)
mtext("Genetic Overlap", side=3, cex=.75)
dev.off()

# Which traits share SNPs:
# Cluster by genetics:
# Map gen to epi loc:
dt = as.dist(1 - genmat)
ht <- hclust(dt, method='ward.D')
cocl <- order.optimal(dt, ht$merge)$order
rn <- names(cocl)[cocl]

epimat = epimat[rn,rn]
tmp = genmat[rn, rn]

diag(tmp) = NA
diag(epimat) = NA
vcol = as.character(ctodf[rn,'COLOR'])
vcol.f = factor(ctodf[rn,'COLOR'])

# Generate breaks:
NCLUST = 20
breaks = calc.breaks(ht, NCLUST, cocl)
acut <- cutree(ht, NCLUST)[cocl]


colramp = col1
# png(paste0(imgpref, 'gwas_genet_cosine_heatmap_bygenet.png'),res=450,units='in',width=18,height=8)
pdf(paste0(imgpref, 'gwas_genet_cosine_heatmap_bygenet.pdf'),width=12,height=6)
layout(matrix(c(1:4), 1, 4),  widths=c(.4,3,.1, 3), heights=3, TRUE)
sp=0.1
par(yaxs='i')
par(xaxs='i')
minor = 1
par(mar=c(sp,rsp,1,0))
image(t(as.matrix(as.numeric(vcol.f))), col=levels(vcol.f), useRaster=T, axes=F)
text(x=par()$usr[1] - 0.01 * diff(par()$usr[1:2]),
     y=seq(0,1,length.out=nrow(epimat)), 
     rownames(epimat), cex=0.1, adj=1, xpd=TRUE)
par(mar=c(sp,sp,1,sp))
# Color sides by the vcol
image(epimat, axes=F, col=c('white', col1), useRaster=TRUE)
mtext("Epigenetic enrichment similarity", side=3, cex=.75)
abline(v=breaks,lty=2,lwd=.25)
abline(h=breaks,lty=2,lwd=.25)
box(lwd=.5)
par(mar=c(sp,0,1,0))
image(t(as.matrix(as.numeric(vcol.f))), col=levels(vcol.f), useRaster=T, axes=F)
par(mar=c(sp,sp,1,sp))
image(1 * (tmp > 0), axes=F, col=c('white', col1), useRaster=TRUE)
abline(v=breaks,lty=2,lwd=.25)
abline(h=breaks,lty=2,lwd=.25)
box(lwd=.5)
mtext("Genetic Overlap", side=3, cex=.75)
dev.off()



# -------------------------------------
# GWAS similarity matrix:
# note: keep all nodes, even empty ones
# -------------------------------------
# Use the normalized node matrix:
Ztmat = redmat %*% ntmat
Ztmat = sweep(Ztmat, 1, apply(Ztmat, 1, sum), '/')

# Try to 1 and jaccard as well?
dt = dist(Ztmat, 'cosine')
dt = dist(Ztmat, 'ejaccard')

dmat = as.matrix(dt)
ht <- hclust(dt, method='ward.D')
cocl <- order.optimal(dt, ht$merge)$order
rn <- names(cocl)[cocl]
mat = dmat[rn,rn]

# Generate breaks:
NCLUST = 20
breaks = calc.breaks(ht, NCLUST, cocl)
acut <- cutree(ht, NCLUST)[cocl]

quant = 0.005
q1 = quantile(mat, quant)
q2 = quantile(mat, 1 - quant)
mat[mat < q1] <- q1
mat[mat > q2] <- q2

colramp=colspec
png(paste0(imgpref, 'gwassim_by_nodes.png'),res=450,units='in',width=20,height=18)
par(mar=c(.5,20,.5,.5))
par(yaxs='i')
par(xaxs='i')
minor = 1
image(mat, axes=F, zlim=c(q1,q2), col=colramp, useRaster=TRUE)  # Raster more efficient?
yt = seq(0,1,length.out=ncol(mat))
xt = par()$usr[3] - 0.001*(par()$usr[4]-par()$usr[3])
yaxlab=rownames(mat)
text(y=yt, x=xt, labels=yaxlab, srt=0, adj=1, xpd=TRUE,cex=.25)
abline(v=breaks,lty=2,lwd=.5)
abline(h=breaks,lty=2,lwd=.5)
dev.off()


