#!/usr/bin/R
# ------------------------------------
# Plot basic motifs enrichment output:
# ------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))
today <- format(Sys.time(), "%m%d%y")
library(ggplot2)
library(viridis)
library(tidyr)

library(rhdf5)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(Matrix)
library(scales)
library(plotrix)
library(PRROC)
library(huxtable)
library(DescTools) # For Gini

imgdir = paste0(img, 'motifs/')
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'raw_')

# df = read.delim(gzfile('motifs/collated.cls.enrich.tsv.gz'), header=T)
df = read.delim(gzfile('motifs/collated.bkgdhs.cls.enrich.tsv.gz'), header=T)
# df = read.delim(gzfile('motifs/collated.bkgdhs.epi.enrich.tsv.gz'), header=T)
map = read.delim('motif_clusters_vierstra2020.tsv', header=T)
colnames(map)[2] = 'motif'
df = merge(df, map[,c('Cluster_ID','motif')])
df$pval_raw[is.na(df$pval_raw)] = 0 
df$pval_ctrl[is.na(df$pval_ctrl)] = 0 
df = df[!is.na(df$pfg),]


# -----------------------------
# Fix the binomial enrichments:
# -----------------------------
calc_binomial_CI = function(ns, n, z, max=T){
    nf = n - ns
    zsq = z^2
    p = (ns + zsq/2) / (n + zsq)
    pm = (z / (n + zsq)) * sqrt((ns * nf) / n + zsq / 4)
    return(p + (2 * max - 1) * pm)}

# Fix such that all enrichment ratios (enrich/deplete) are more conservative:
df$rc = log2((df$pfg / df$pbg) / (df$cfg / df$cbg))
pind = which(df$rc >= 0)
nind = which(df$rc < 0)
df$pfrac[pind] = sapply(pind, function(i){calc_binomial_CI(df$pfg[i], df$pbg[i], 1.5, max=F)})
df$cfrac[pind] = sapply(pind, function(i){calc_binomial_CI(df$cfg[i], df$cbg[i], 1.5, max=T)})
df$pfrac[nind] = sapply(nind, function(i){calc_binomial_CI(df$pfg[i], df$pbg[i], 1.5, max=T)})
df$cfrac[nind] = sapply(nind, function(i){calc_binomial_CI(df$cfg[i], df$cbg[i], 1.5, max=F)})
df$pfrac[is.na(df$pfrac)] = 0
df$cfrac[is.na(df$cfrac)] = 0
# log2 fold enrichment:
df$ratio_ctrl = log2(df$pfrac / df$cfrac)

# Re-do the pval to acct. for depletion. ** (two-sided pval at 0.05)
enrp = apply(df[,c('pfg','pbg','cfg','cbg')], 1, 
             function(y){phyper(q=y[1]-1, m=y[3], n=y[4] - y[3], k=y[2], lower.tail=FALSE)})
depp = apply(df[,c('pfg','pbg','cfg','cbg')], 1, 
             function(y){phyper(q=y[1], m=y[3], n=y[4] - y[3], k=y[2], lower.tail=TRUE)})
df$pval_ctrl = apply(cbind(-log10(enrp), -log10(depp)),1, max)

# For raw:
ntot = sum(unique(df[,c('set','regtot')])$regtot)
df$ratio_raw = log2((df$pfg / df$pbg) / (df$regtot / ntot))


# Keep motifs with p < 0.01 and log2FC > 1.5
cutoff = 1.5
pcutoff = -log10(0.005)
sdf = df[df$pval_ctrl >= pcutoff & abs(df$ratio_ctrl) >= 1.5,]
rdf = df[df$pval_raw > 2 & abs(df$ratio_raw) > 1.5,]
top.ctrl = aggregate(ratio_ctrl ~ Cluster_ID, sdf, function(x){x[which.max(abs(x))]}) 
top.raw = aggregate(ratio_raw ~ Cluster_ID, rdf, function(x){x[which.max(abs(x))]}) 
cdf = aggregate(motif ~ Cluster_ID, merge(sdf, top.ctrl), function(x){x[1]})
rdf = aggregate(motif ~ Cluster_ID, merge(rdf, top.raw), function(x){x[1]})
ctrl.motifs = cdf$motif
raw.motifs = rdf$motif

# --------------------------------------
# Plot the controlled motif enrichments:
# --------------------------------------
rwide = spread(df[df$motif %in% ctrl.motifs,
               c('set','motif','ratio_ctrl')], set, ratio_ctrl, fill=0)
pwide = spread(df[df$motif %in% ctrl.motifs,
               c('set','motif','pval_ctrl')], set, pval_ctrl, fill=0)
rmat = as.matrix(rwide[,-1])
pmat = as.matrix(pwide[,-1])
rownames(rmat) = rwide[,1]
rownames(pmat) = pwide[,1]
pmat[pmat < 2] = 0
pmat[pmat >= 2] = 1
rmat[is.na(rmat)] = 0
# zlim = 2
zlim = 1.5
rmat[rmat < -zlim] = -zlim
rmat[rmat > zlim] = zlim
rmat = rmat * pmat[rownames(rmat), colnames(rmat)]
print(dim(rmat))
rmat = rmat[, apply(rmat,2, max) >= 1.5]

# Cluster matrix:
dt <- dist(rmat, 'euclidean')
roword <- hclust(dt)$order
dt <- dist(t(rmat), 'euclidean')
colord <- hclust(dt)$order
rmat = rmat[roword, colord]
rmat = diag.mat2(rmat)[[1]]
reord = t(rmat)

png(paste0('test_ctrl_heatmap.png'), res=400, width=7, height=7.0 * ncol(reord) / nrow(reord), units='in')
par(mar=c(.25,4,.25,.25))
image(reord, axes=FALSE,col=rev(colrb), zlim=c(-zlim, zlim), useRaster=T)
# grid(nx=nrow(reord), ny=ncol(reord),col='grey',lty='solid',lwd=.25)
text(y=seq(0,1,length.out=ncol(reord)), 
     x=par()$usr[1]-0.003*(par()$usr[2]-par()$usr[1]),
     labels=colnames(reord), srt=0, adj=1, xpd=TRUE,cex=.15)
text(x=seq(0,1,length.out=nrow(reord)), 
     y=par()$usr[3]-0.002*(par()$usr[4]-par()$usr[3]),
     labels=rownames(reord), srt=90, adj=1, xpd=TRUE,cex=.15)
box(lwd=.25)
dev.off()

subdf = df[df$motif %in% ctrl.motifs,]
subdf = subdf[subdf$pval_ctrl >= pcutoff,]

# ------------------------------------------
# Plot the ALL controlled motif enrichments:
# ------------------------------------------
rwide = spread(df[, c('set','motif','ratio_ctrl')], set, ratio_ctrl, fill=0)
pwide = spread(df[, c('set','motif','pval_ctrl')], set, pval_ctrl, fill=0)
rmat = as.matrix(rwide[,-1])
pmat = as.matrix(pwide[,-1])
rownames(rmat) = rwide[,1]
rownames(pmat) = pwide[,1]
pmat[pmat < pcutoff] = 0
pmat[pmat >= pcutoff] = 1
rmat[is.na(rmat)] = 0
# zlim = 2
zlim = 1.5
rmat[rmat < -zlim] = -zlim
rmat[rmat > zlim] = zlim
rmat = rmat * pmat[rownames(rmat), colnames(rmat)]
print(dim(rmat))
rmat = rmat[, apply(rmat,2, max) >= 1.5]
rmat = rmat[apply(rmat,1, max) >= 1.5,]

# Cluster matrix:
dt <- dist(rmat, 'euclidean')
roword <- hclust(dt)$order
dt <- dist(t(rmat), 'euclidean')
colord <- hclust(dt)$order
rmat = rmat[roword, colord]
rmat = diag.mat2(rmat)[[1]]
reord = t(rmat)

png(paste0('test_ctrl_all_heatmap.png'), res=400, width=7, height=7.0 * ncol(reord) / nrow(reord), units='in')
par(mar=c(.25,4,.25,.25))
image(reord, axes=FALSE,col=rev(colrb), zlim=c(-zlim, zlim), useRaster=T)
# grid(nx=nrow(reord), ny=ncol(reord),col='grey',lty='solid',lwd=.25)
text(y=seq(0,1,length.out=ncol(reord)), 
     x=par()$usr[1]-0.003*(par()$usr[2]-par()$usr[1]),
     labels=colnames(reord), srt=0, adj=1, xpd=TRUE,cex=.15)
text(x=seq(0,1,length.out=nrow(reord)), 
     y=par()$usr[3]-0.002*(par()$usr[4]-par()$usr[3]),
     labels=rownames(reord), srt=90, adj=1, xpd=TRUE,cex=.15)
box(lwd=.25)
dev.off()


# -------------------------------
# Plot the raw motif enrichments:
# -------------------------------
rwide = spread(df[df$motif %in% raw.motifs,
               c('set','motif','ratio_raw')], set, ratio_raw, fill=0)
rmat = as.matrix(rwide[,-1])
rownames(rmat) = rwide[,1]
rmat[is.na(rmat)] = 0
zlim = 2
rmat[rmat < -zlim] = -zlim
rmat[rmat > zlim] = zlim

# Cluster matrix:
dt <- dist(rmat, 'euclidean')
roword <- hclust(dt)$order
dt <- dist(t(rmat), 'euclidean')
colord <- hclust(dt)$order
reord = rmat[roword, colord]
reord = t(reord)

png(paste0('test_raw_heatmap.png'), res=400, width=7, height=7.0 * ncol(reord) / nrow(reord), units='in')
par(mar=c(1,4,.25,.25))
image(reord, axes=FALSE,col=rev(colrb), zlim=c(-zlim, zlim), useRaster=TRUE)
# grid(nx=nrow(reord), ny=ncol(reord),col='grey',lty='solid',lwd=.25)
text(y=seq(0,1,length.out=ncol(reord)), 
     x=par()$usr[1]-0.003*(par()$usr[2]-par()$usr[1]),
     labels=colnames(reord), srt=0, adj=1, xpd=TRUE,cex=.15)
text(x=seq(0,1,length.out=nrow(reord)), 
     y=par()$usr[3]-0.002*(par()$usr[4]-par()$usr[3]),
     labels=rownames(reord), srt=90, adj=1, xpd=TRUE,cex=.15)
box(lwd=.25)
dev.off()


# ----------------------------------------------------------
# Get enrichments over the full enhancer background instead:
# ----------------------------------------------------------
totclsdf = aggregate(cbind(pfg, cfg) ~ motif, df, sum)
names(totclsdf) = c('motif','pbgdhs','cbgdhs')
df = merge(df, totclsdf, all.x=TRUE)

totregdf = aggregate(regtot ~ motif, df, sum)
totdhs = totregdf$regtot[1]

# Recalc ratios:
# Fix such that all enrichment ratios (enrich/deplete) are more conservative:
df$rc = log2((df$pfg / df$pbgdhs) / (df$cfg / df$cbgdhs))
pind = which(df$rc >= 0)
nind = which(df$rc < 0)
df$pfrac[pind] = sapply(pind, function(i){calc_binomial_CI(df$pfg[i], df$pbgdhs[i], 1.5, max=F)})
df$cfrac[pind] = sapply(pind, function(i){calc_binomial_CI(df$cfg[i], df$cbgdhs[i], 1.5, max=T)})
df$pfrac[nind] = sapply(nind, function(i){calc_binomial_CI(df$pfg[i], df$pbgdhs[i], 1.5, max=T)})
df$cfrac[nind] = sapply(nind, function(i){calc_binomial_CI(df$cfg[i], df$cbgdhs[i], 1.5, max=F)})
df$pfrac[is.na(df$pfrac)] = 0
df$cfrac[is.na(df$cfrac)] = 0
# log2 fold enrichment:
df$ratio_ctrl = log2(df$pfrac / df$cfrac)

# TODO: Recalculate the ctrl pvalue:

sdf = df[df$pval_ctrl > 2 & abs(df$ratio_ctrl) > 1.5,]
top.ctrl = aggregate(ratio_ctrl ~ Cluster_ID, sdf, function(x){x[which.max(abs(x))]}) 
cdf = aggregate(motif ~ Cluster_ID, merge(sdf, top.ctrl), function(x){x[1]})
ctrl.motifs = cdf$motif

# --------------------------------------
# Plot the controlled motif enrichments:
# --------------------------------------
rwide = spread(df[df$motif %in% ctrl.motifs,
               c('set','motif','ratio_ctrl')], set, ratio_ctrl, fill=0)
rmat = as.matrix(rwide[,-1])
rownames(rmat) = rwide[,1]
rmat[is.na(rmat)] = 0
# zlim = 2
zlim = 1.5
rmat[rmat < -zlim] = -zlim
rmat[rmat > zlim] = zlim
print(dim(rmat))
rmat = rmat[, apply(rmat,2, max) >= 1.5]

# Cluster matrix:
# dt <- dist(rmat, 'euclidean')
# roword <- hclust(dt)$order
# dt <- dist(t(rmat), 'euclidean')
# colord <- hclust(dt)$order
# rmat = rmat[roword, colord]
rmat = diag.mat2(rmat)[[1]]
reord = t(rmat)

png(paste0('test_ctrl_heatmap2.png'), res=400, width=7, height=7.0 * ncol(reord) / nrow(reord), units='in')
par(mar=c(.25,4,.25,.25))
image(reord, axes=FALSE,col=rev(colrb), zlim=c(-zlim, zlim), useRaster=T)
# grid(nx=nrow(reord), ny=ncol(reord),col='grey',lty='solid',lwd=.25)
text(y=seq(0,1,length.out=ncol(reord)), 
     x=par()$usr[1]-0.003*(par()$usr[2]-par()$usr[1]),
     labels=colnames(reord), srt=0, adj=1, xpd=TRUE,cex=.15)
text(x=seq(0,1,length.out=nrow(reord)), 
     y=par()$usr[3]-0.002*(par()$usr[4]-par()$usr[3]),
     labels=rownames(reord), srt=90, adj=1, xpd=TRUE,cex=.15)
box(lwd=.25)
dev.off()


# ----------------------
# Recalc ratios for raw:
# ----------------------
# Fix such that all enrichment ratios (enrich/deplete) are more conservative:
df$rcraw = log2((df$pfg / df$pbgdhs) / (df$regtot / totdhs))
pind = which(df$rcraw >= 0)
nind = which(df$rcraw < 0)
df$pfracraw = df$pfrac
df$cfracraw = df$cfrac
df$pfracraw[pind] = sapply(pind, function(i){calc_binomial_CI(df$pfg[i], df$pbgdhs[i], 1.5, max=F)})
df$cfracraw[pind] = sapply(pind, function(i){calc_binomial_CI(df$regtot[i], totdhs, 1.5, max=T)})
df$pfracraw[nind] = sapply(nind, function(i){calc_binomial_CI(df$pfg[i], df$pbgdhs[i], 1.5, max=T)})
df$cfracraw[nind] = sapply(nind, function(i){calc_binomial_CI(df$regtot[i], totdhs[i], 1.5, max=F)})
df$cfracraw = df$regtot / totdhs
df$pfracraw[is.na(df$pfracraw)] = 0
df$cfracraw[is.na(df$cfracraw)] = 0
# log2 fold enrichment:
df$ratio_raw = log2(df$pfracraw / df$cfracraw)

# TODO: Recalculate the ctrl pvalue:

# Keep motifs with log2FC > 1.5
rdf = df[abs(df$ratio_raw) > 1.5,]
top.raw = aggregate(ratio_raw ~ Cluster_ID, rdf, function(x){x[which.max(abs(x))]}) 
rdf = aggregate(motif ~ Cluster_ID, merge(rdf, top.raw), function(x){x[1]})
raw.motifs = rdf$motif


# -------------------------------
# Plot the raw motif enrichments:
# -------------------------------
rwide = spread(df[df$motif %in% raw.motifs,
               c('set','motif','ratio_raw')], set, ratio_raw, fill=0)
rmat = as.matrix(rwide[,-1])
rownames(rmat) = rwide[,1]
rmat[is.na(rmat)] = 0
zlim = 2
rmat[rmat < -zlim] = -zlim
rmat[rmat > zlim] = zlim

# Cluster matrix:
dt <- dist(rmat, 'euclidean')
roword <- hclust(dt)$order
dt <- dist(t(rmat), 'euclidean')
colord <- hclust(dt)$order
rmat = rmat[roword, colord]
rmat = diag.mat2(rmat)[[1]]
reord = t(rmat)

png(paste0('test_raw_heatmap2.png'), res=400, width=7, height=7.0 * ncol(reord) / nrow(reord), units='in')
par(mar=c(1,4,.25,.25))
image(reord, axes=FALSE,col=rev(colrb), zlim=c(-zlim, zlim), useRaster=TRUE)
# grid(nx=nrow(reord), ny=ncol(reord),col='grey',lty='solid',lwd=.25)
text(y=seq(0,1,length.out=ncol(reord)), 
     x=par()$usr[1]-0.003*(par()$usr[2]-par()$usr[1]),
     labels=colnames(reord), srt=0, adj=1, xpd=TRUE,cex=.15)
text(x=seq(0,1,length.out=nrow(reord)), 
     y=par()$usr[3]-0.002*(par()$usr[4]-par()$usr[3]),
     labels=rownames(reord), srt=90, adj=1, xpd=TRUE,cex=.15)
box(lwd=.25)
dev.off()


