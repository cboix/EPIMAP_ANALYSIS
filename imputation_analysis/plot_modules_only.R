#!/usr/bin/R
# ------------------------------------------
# Plot modules from CHMM/other clustering
# 1. Plot the centers matrix
# 2. Plot the margins 
# (Number of epigenomes per cluster, vice versa)
# 3. Plot and order/cluster the adjacency matrix
#       for the epigenomes
# ------------------------------------------
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
library(dplyr)

# Arguments:
filepref = 'cls_merge2_wH3K27ac100_300'

tagline = 'ChromHMM Enhancers'
imgdir = paste0(img, "clusters/") 
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need clusters filename")
} else {        
    filepref = args[1]
    tagline = args[2]
    if (length(args) > 2){ 
        imgdir = args[3] 
    } else { 
        imgdir = paste0(img, "clusters/") 
    }
}

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){
    c(filepref, tagline, imgdir) }
source(paste0(bindir, 'load_modules_data.R'))

# Prefix:
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'clusters_')

parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}

# ---------------------
# Plot by frozen order:
# ---------------------
png(paste0(imgpref,'centers_orddist.png'),res=300,units='in',width=14,height=12)
dist.order = plot_clusters(rn, acut.nam, lablist.orddist)
dev.off()

# -------------------------------
# Plot ordered by centers matrix:
# -------------------------------
png(paste0(imgpref,'centers_matdist.png'),res=300,units='in',width=14,height=12)
dist.order = plot_clusters(rnmat, acutmat.nam, lablist.mat)
dev.off()

# -------------------------------------------
# Plot ordered by groups, then centers matrix
# -------------------------------------------
png(paste0(imgpref,'centers_groupord.png'),res=300,units='in',width=14,height=12)
dist.order = plot_clusters(rngroup, acutgroup.nam, lablist.group, plot.rect=TRUE)
dev.off()

# -------------------------------------------
# Plot ordered by groups, then centers matrix
# -------------------------------------------
png(paste0(imgpref,'centers_groupord.png'),res=300,units='in',width=14,height=12)
dist.order = plot_clusters(rngroup, acutgroup.nam, lablist.group, plot.rect=TRUE)
dev.off()

png('~/test.png',res=300,units='in',width=8,height=7)
png(paste0(imgpref,'centers_groupord_medium.png'),res=300,units='in',width=8,height=7)
dist.order = plot_clusters_vsmall(rngroup, acutgroup.nam, lablist.group, plot.rect=TRUE)
dev.off()

png(paste0(imgpref,'centers_groupord_small.png'),res=300,units='in',width=7.75,height=5)
dist.order = plot_clusters_vsmall(rngroup, acutgroup.nam, lablist.group, plot.rect=TRUE, small=TRUE, count.barplot=FALSE)
dev.off()

pdf(paste0(imgpref,'centers_groupord_small.pdf'), width=7.75, height=5)
dist.order = plot_clusters_vsmall(rngroup, acutgroup.nam, lablist.group, plot.rect=TRUE, small=TRUE, count.barplot=FALSE)
dev.off()


png(paste0(imgpref,'centers_groupord_tiny.png'),res=300,units='in',width=7.75,height=4.25)
dist.order = plot_clusters_vsmall(rngroup, acutgroup.nam, lablist.group, plot.rect=TRUE, tiny=TRUE, count.barplot=FALSE)
dev.off()

pdf(paste0(imgpref,'centers_groupord_tiny.pdf'), width=7.75, height=4.25)
dist.order = plot_clusters_vsmall(rngroup, acutgroup.nam, lablist.group, plot.rect=TRUE, tiny=TRUE, count.barplot=FALSE)
dev.off()


# -----------------------
# Add the three palettes:
# -----------------------
# 1. Cluster centers
# 2. GO terms
# 3. Motifs
# -----------------------
# Update palette and legend for pvalues:
cls.col_fun = function(x, pal=rev(col1)){
    palette = rev(pal)
    bin <- cut(x, seq(0, 100, length.out=length(palette)), include.lowest=T) 
    palette[bin] }
cls.legend = Legend(at = seq(0, 100, 25), 
                    labels_gp = gpar(fontsize=5),
                    title_gp = gpar(fontsize=5, fontface='bold'),
                    col_fun=cls.col_fun, title_position = "topleft", 
                    title='Inclusion %', direction = 'vertical')

go.col_fun = function(x){
    x[x < 2] = 1
    x[x < 3 & x >= 2] = 2
    x[x < 4 & x >= 3] = 3
    x[x >= 4] = 4
    rcols = c('white','yellow','orange','red')
    rcols[x] }
go.legend = Legend(at = 1:5,
                   labels = c('0','2','3','4','>4'),
                   labels_gp = gpar(fontsize=5),
                   title_gp = gpar(fontsize=5, fontface='bold'),
                   col_fun=go.col_fun, title_position = "topleft", 
                   title='-log10p', direction = 'vertical')

mot.col_fun = function(x, pal=colrb){
    palette = rev(pal)
    bin <- cut(x, seq(-1.5, 1.5, length.out=length(palette)), include.lowest=T) 
    palette[bin] }
mot.legend = Legend(at = seq(-1.5, 1.5, .5), 
                    labels_gp = gpar(fontsize=5),
                    title_gp = gpar(fontsize=5, fontface='bold'),
                    col_fun=mot.col_fun, title_position = "topleft", 
                    title='log2\nEnrich.', direction = 'vertical')

plegend = packLegend(cls.legend, go.legend, mot.legend)

pdf(paste0(imgpref,'centers_groupord_tiny_wlegend.pdf'), width=7.75, height=4.25)
dist.order = plot_clusters_vsmall(rngroup, acutgroup.nam, lablist.group, plot.rect=TRUE, tiny=TRUE, count.barplot=FALSE)
draw(plegend, x = unit(0.25,'in'), y=unit(4,'in'), just = "top")
dev.off()



# ----------------------------------
# See specificity of module centers:
# ----------------------------------
Z = 1 * (centers > 0.25)
rn = rownames(Z)
labels = meta[rn, 'GROUP']
marg = apply(Z, 2, sum)

zdf = data.frame(Z)
zdf$GROUP = labels
zw = gather(zdf, cls, val, - GROUP)
zw = zw[zw$val > 0,]
mdf = aggregate(val ~ cls, zw, sum)
names(mdf)[2] = 'total'

cdf = aggregate(val ~ cls + GROUP, zw, sum)
gdf = aggregate(GROUP ~ cls, cdf, length)
cdf = merge(cdf, mdf)
cdf$ratio = cdf$val / cdf$total
cdf = aggregate(ratio ~ cls, cdf, max)
sum(cdf$ratio >= 0.9)
sum(cdf$ratio >= 0.5)

# Numbers for manuscript
ind = which(marg >= 150)
rind = which(marg < 150)
mean(marg[rind]) / 833 * 100
mean(marg[ind]) / 833 * 100

range(gdf$GROUP)
sum(gdf$GROUP >= .1 * 33)
mean(gdf$GROUP[gdf$cls %in% names(ind)]) / 33
mean(mdf$total[mdf$cls %in% names(ind)]) / 833

sum(counts[[tagline]][ind])

# Numbers of enriched for adult/embryonic (that aren't single-sample):
a.ls = names(which(enrichmat[,'lifestage'] == 'adult'))
e.ls = names(which(enrichmat[,'lifestage'] == 'embryonic'))
length(a.ls)
length(e.ls)

a.mls = names(which(metamat[,'lifestage'] == 'adult'))
e.mls = names(which(metamat[,'lifestage'] == 'embryonic'))

# Avg. size:
a.tot = mdf$total[mdf$cls %in% a.ls[!(a.ls %in% names(ind))]]
e.tot = mdf$total[mdf$cls %in% e.ls[!(e.ls %in% names(ind))]]
mean(a.tot)
mean(e.tot)




# -------------------------------------------
# Plot hypergeometric enrichments of centers:
# -------------------------------------------
# FLIP: 
hgmat = hgmat[rev(1:nrow(hgmat)),]
clsord = dist.order[[1]]
vbreaks = dist.order[[2]]
sublab = sapply(rownames(hgmat), function(x) sub(".*\n","",x))
subgroup = sapply(rownames(hgmat), function(x) sub("\n.*","",x))
subgn = as.numeric(factor(subgroup, levels=unique(subgroup)))
mot.col_fun = function(x, palette=col1){
    bin <- cut(x, seq(0, 5, length.out=length(palette)), include.lowest=T) 
    palette[bin] }
mot.legend = Legend(at=seq(0,5, 1), 
                    labels_gp = gpar(fontsize=5),
                    title_gp = gpar(fontsize=5, fontface='bold'),
                    col_fun=mot.col_fun, title_position = "topleft", 
                    legend_height=unit(.5,'in'),
                    title='-log10p', direction = 'vertical')
plegend = packLegend(mot.legend)

# NOTE: orddist is not alpha ordered, may not match group lab order
# pdf(paste0(imgpref, 'hg_enrich_clusters.pdf'), width=8, height=3.5)
pdf(paste0(imgpref, 'hg_enrich_clusters.pdf'), width=8, height=3)
# png(paste0(imgpref, 'hg_enrich_clusters.png'), res=250, width=20, height=6, units='in')
sp = 0.1
par(mar=c(1.25,4,sp,sp))
image(t(hgmat[,clsord]), col=col, axes=F, useRaster=T)
text(x=seq(0,1,length.out=ncol(hgmat)), 
     y=par()$usr[3]-0.005*(par()$usr[4]-par()$usr[3]),
     labels=clsord,  # colnames(hgmat), 
     srt=90, adj=1, xpd=TRUE,cex=.2)
text(y=seq(0,1,length.out=nrow(hgmat)), 
     x=par()$usr[1]-0.002*(par()$usr[2]-par()$usr[1]), 
     labels=sublab, srt=0, adj=1, xpd=TRUE,cex=.3, col=ifelse(subgn %% 2 == 1, 'black','grey50'))
# TODO add colors?
cols=(enrbreaks -.5) / (nrow(hgmat) - 1)
abline(h=cols,col='black',lty='dashed',lwd=.25)
rows=(seq(0, ncol(hgmat), 3) -.5) / (ncol(hgmat) - 1)
abline(v=vbreaks,col='black',lty='dashed',lwd=.25)
box(lwd=.5)
mtext(paste0('Covariate Levels'), side=2, cex=.5, line=3)
mtext(paste0(tagline, ' Modules'), side=1, cex=.5, line=0.25)
draw(plegend, x = unit(.20,'in'), y=unit(2.95,'in'), just = "top")
dev.off()


# -----------------------------------
# Plot the margins of cluster matrix:
# -----------------------------------
zmat = cmat
zmat[zmat > threshold] = 1
zmat[zmat < threshold] = 0
zmat2 = cmat
zmat2[zmat2 > .25] = 1
zmat2[zmat2 < .25] = 0
# Margins:
ncls = apply(zmat, 1, sum)
nepi = apply(zmat, 2, sum)
ncls2 = apply(zmat2, 1, sum)
nepi2 = apply(zmat2, 2, sum)

# png(paste0(imgpref, 'margins.png'), res=250, width=8, height=4, units='in')
pdf(paste0(imgpref, 'margins.pdf'), width=8, height=4)
layout(matrix(c(1:4),2,2), TRUE)
par(yaxs="i")
par(xaxs="i")
par(mar=c(3.5,3.5,2,.5))
hist(ncls, 50, col='darkgrey', border='white', xlim=c(0,max(ncls)), 
     ylab='', xlab='', main='')
# box()
mtext('# Epigenomes', side=2,line=2)
mtext('# Clusters', side=1,line=2)
mtext('Centers Thresholded at 10%', side=3,line=.25)
par(mar=c(3.5,3.5,0,.5))
hist(nepi, 50, col='darkgrey', border='white', xlim=c(0,max(nepi)),
     ylab='', xlab='', main='')
# box()
mtext('# Epigenomes', side=1,line=2)
mtext('# Clusters', side=2,line=2)
par(mar=c(3.5,3.5,2,.5))
hist(ncls2, 50, col='darkgrey', border='white',  xlim=c(0,max(ncls2)), 
     ylab='', xlab='', main='')
# box()
mtext('# Epigenomes', side=2,line=2)
mtext('# Clusters', side=1,line=2)
mtext('Centers Thresholded at 25%', side=3,line=.25)
par(mar=c(3.5,3.5,0,.5))
hist(nepi2, 50, col='darkgrey', border='white',  xlim=c(0,max(nepi2)), 
     ylab='', xlab='', main='')
mtext('# Epigenomes', side=1,line=2)
mtext('# Clusters', side=2,line=2)
# box()
dev.off()



# ------------------------------------------
# Plot the ejaccard dist between epigenomes:
# ------------------------------------------
method = 'jaccard'
zmat = cmat
zmat[zmat > threshold] = 1
zmat[zmat < threshold] = 0
dt = dist(zmat, method)
# dt = dist(cmat, method)
dmat = as.matrix(dt)
ht <- hclust(dt, method='ward.D')
cocl <- order.optimal(dt, ht$merge)$order
rnmat  <- names(cocl)[cocl]

# Generate breaks:
NCLUST = 20
breaks = calc.breaks(ht, NCLUST, cocl)
acutmat <- cutree(ht, NCLUST)[cocl]
dist.breaks = calc.breaks.acut(acutmat)
acutmat.nam = acutmat
names(acutmat.nam) = rnmat
reord = dmat[rnmat, rnmat]
palette = rev(viridis(100))

# Labels in figure:
labelsmat = meta[rnmat, 'GROUP']
faclabels = as.matrix(as.numeric(labelsmat))
colset = as.character(rdcol$COLOR)
fix = rnmat[is.na(labelsmat)]
print(length(fix))
lablist = label.runs(faclabels, labelsmat, rdcol)
palette = colryb

# Plot matrix:
# png(paste0(imgpref, 'matdist.png'), res=250, width=12, height=10, units='in')
pdf(paste0(imgpref, 'matdist.pdf'), width=12, height=10)
layout(matrix(c(1:4),2,2), heights=c(8,.5), widths=c(1.5, 8), TRUE)
par(yaxs="i")
par(xaxs="i")
par(mar=c(0.25, 6, 2, 0))
meta.image(metamat[rnmat,5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
text(y=lablist[[1]],
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lablist[[2]], srt=0, adj=1, xpd=TRUE,cex=.8, col=lablist[[3]])
abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
par(mar=c(0.25, 6, 0, 0))
image(t(as.matrix(metamat) == ""), col='white', axes=F, useRaster=TRUE)
metaclass = sapply(rev(colnames(metamat)), capitalize)
text(x=seq(0,1, length.out=ncol(metamat)),
     y=par()$usr[4]-0.01*(par()$usr[4]-par()$usr[3]), 
     labels=metaclass, srt=90, adj=1, xpd=TRUE, cex=.7)
par(mar=c(0.25, .25, 2, 0))
plot.cov(reord, clamp=TRUE, palette=palette, breaks=breaks)
mtext(paste0("Modules shared (Jaccard similarity) in ", tagline), side=3, line=.25, cex=1.3)
dev.off()


# -------------------------------------------------
# Look at the enrichments of the cancer epigenomes:
# -------------------------------------------------
# Ask - samples with which groups are co-enriched in modules with cancer samples?
# cent.adj = centers %*% t(centers)

# cent.adj = (centers > .1) %*% t((centers > .1))
cent.adj = (centers[,rind]) %*% t(centers[,rind])
rn.grp = as.character(meta[rownames(cent.adj),'GROUP'])
uqgrp = odf$GROUP
uqgrp = uqgrp[!(uqgrp %in% c('Multiple', 'Cancer'))]
# tform = make.tform(rn.grp, u=uqgrp, norm=TRUE)
tform = make.tform(rn.grp, u=uqgrp, norm=TRUE)
cent.grp = cent.adj %*% tform
cind = which(rn.grp == 'Cancer')

ca.grp = cent.grp[cind,]
ca.grp = sweep(ca.grp, 1, rowSums(ca.grp),'/')
# diag.mat2(ca.grp)
# image(scale(t(ca.grp)))

ca.max = colnames(ca.grp)[apply(ca.grp, 1, which.max)]
ca.2nd = as.character(meta[rownames(ca.grp),'SECONDARY'])

# Reorder:
mat = ca.grp
dt = dist(mat)
ht <- hclust(dt, method='ward.D')
roword <- order.optimal(dt, ht$merge)$order
dt = dist(t(mat))
ht <- hclust(dt, method='ward.D')
colord <- order.optimal(dt, ht$merge)$order
plt.mat = ca.grp[roword,colord]
rn = rownames(plt.mat)
cn = colnames(plt.mat)

scale = 15
h = (dim(ca.grp)[1] + 7) / scale
w = (dim(ca.grp)[2] + 11) / scale

pdf(paste0(imgpref,'cancer_assignment.pdf'), width=w, height=h)
# png(paste0(imgpref,'cancer_assignment.png'), res=450, units='in', width=w, height=h)
# par(mar=c(4,4,0,0))
par(mar=c(4,4,0,0))
image(t(plt.mat), axes=F, useRaster=T)
box(lwd=.5)
xat = seq(0,1, length.out=ncol(ca.grp))
yat = seq(0,1, length.out=nrow(ca.grp))
text(y=par()$usr[3] - .005 * diff(par()$usr[3:4]), 
     x=xat, labels=colnames(plt.mat), xpd=TRUE, srt=90, adj=1, cex=.4)
text(x=par()$usr[1] - .005 * diff(par()$usr[1:2]), 
     y=yat, labels=meta[rownames(plt.mat), 'infoline'], xpd=TRUE, adj=1, cex=.25)
# HIGHLIGHT WHERE THE TOP ENRICHMENT SHOULD BE.
gind = sapply(meta[rownames(plt.mat), 'SECONDARY'], function(x){which(colnames(plt.mat) == x)})
points(xat[gind], yat, pch=22, cex=1, lwd=.5)
# Max:
mind = apply(plt.mat, 1, which.max)
points(xat[mind], yat, pch='*', cex=1, lwd=.25)
dev.off()



scale = 20
widths = c(11, 1, dim(ca.grp)[2])
heights = c(dim(ca.grp)[1], 1, 7) 
h = sum(heights) / scale
w = sum(widths) / scale 

pdf(paste0(imgpref,'cancer_assignment_withcolors.pdf'), width=w, height=h)
# png(paste0(imgpref,'cancer_assignment_withcolors.png'), res=450, units='in', width=w, height=h)
layout(matrix(1:9, nrow=3,ncol=3, byrow=TRUE), widths=widths, heights=heights, TRUE)
xat = seq(0,1, length.out=ncol(ca.grp))
yat = seq(0,1, length.out=nrow(ca.grp))
par(mar=rep(0,4))
image(t(plt.mat), col='white', axes=F, useRaster=T)
text(x=par()$usr[2] - .001 * diff(par()$usr[1:2]), 
     y=yat, labels=meta[rownames(plt.mat), 'infoline'], xpd=TRUE, adj=1, cex=.25)
# COLORS:
plt.2nd = meta[rownames(plt.mat), 'SECONDARY']
f2nd = factor(plt.2nd)
image(t(as.numeric(f2nd)), col=colvals$group[levels(f2nd)], axes=F, useRaster=T)
# par(mar=c(4,4,0,0))
# image(t(plt.mat), axes=F, useRaster=T)
image(t(plt.mat), axes=F, useRaster=T)
box(lwd=.5)
# HIGHLIGHT WHERE THE TOP ENRICHMENT SHOULD BE.
gind = sapply(plt.2nd, function(x){which(colnames(plt.mat) == x)})
points(xat[gind], yat, pch=22, cex=1, lwd=.5)
# Max:
mind = apply(plt.mat, 1, which.max)
points(xat[mind], yat, pch='*', cex=1, lwd=.25)
plot(1,1, type='n', axes=F)
plot(1,1, type='n', axes=F)
frow = factor(colnames(plt.mat))
image(t(t(as.numeric(frow))), col=colvals$group[levels(frow)], axes=F, useRaster=T)
plot(1,1, type='n', axes=F)
plot(1,1, type='n', axes=F)
frow = factor(colnames(plt.mat))
image(t(t(as.numeric(frow))), col='white', axes=F, useRaster=T)
text(y=par()$usr[4] - .005 * diff(par()$usr[3:4]), 
     x=xat, labels=colnames(plt.mat), xpd=TRUE, srt=90, adj=1, cex=.4)
dev.off()


# Add legend:
col_fun = function(x, pal=hcl.colors(12, "YlOrRd")){
    palette = rev(pal)
    bin <- cut(x, seq(0, max(plt.mat), length.out=length(palette)), include.lowest=T) 
    palette[bin] }

clegend = Legend(at = seq(0, .15, .05), 
                    labels_gp = gpar(fontsize=5),
                    title_gp = gpar(fontsize=5, fontface='bold'),
                    col_fun=col_fun, title_position = "topleft", 
                    title='Norm.\nAdj.', direction = 'vertical')


pdf(paste0(imgpref,'cancer_assign_legend.pdf'), width=4, height=4)
draw(clegend, x = unit(0.25,'in'), y=unit(2,'in'), just = "top")
dev.off()


# ------------------------------------------------------------
# Find closest non-cancer sample in the raw imputed distances:
# ------------------------------------------------------------
# Load specific distance matrices:
fixedstate = TRUE
nregions = 0 
dataset = 'spearman'

if (fixedstate){
    commandArgs <- function(trailingOnly=TRUE){ c(nregions, dataset) }
    source(paste0(bindir, 'load_region_distance_matrices.R'))
    setprefix = paste('region',nregions, dataset,'distances_', sep="_")
} else { 
    source(paste0(bindir, 'load_distance_matrices.R'))
    setprefix = 'distances_'
    metric = 'genome-wide correlation'
}


# Which is the closest sample:
fna = full
diag(fna) = 1
rn.full = as.character(meta[rownames(fna),'GROUP'])
cind = which(rn.full == 'Cancer')
ncind = which(rn.full != 'Cancer')
cmin = apply(fna[cind,ncind], 1, which.min)
min.samp = rownames(fna)[ncind[cmin]]
ca.samp = rownames(fna)[cind]

meta[min.samp,'GROUP'] == meta[ca.samp,'SECONDARY']
cbind(as.character(meta[min.samp,'GROUP']), meta[ca.samp,'SECONDARY'])


# Average on groups:
rn.full = as.character(meta[rownames(full),'GROUP'])
tform = make.tform(rn.full, u=uqgrp, norm=TRUE)
full.grp = full %*% tform
cind = which(rn.full == 'Cancer')
ca.grp = full.grp[cind,]

# Reorder:
rn2 = rn
rn2 = rn2[rn2 %in% rownames(ca.grp)]
mat = ca.grp
# dt = dist(mat)
# ht <- hclust(dt, method='ward.D')
# roword <- order.optimal(dt, ht$merge)$order
# dt = dist(t(mat))
# ht <- hclust(dt, method='ward.D')
# colord <- order.optimal(dt, ht$merge)$order
plt.mat = ca.grp[rn2,cn]

scale = 15
h = (dim(ca.grp)[1] + 7) / scale
w = (dim(ca.grp)[2] + 11) / scale

pdf(paste0(imgpref,'cancer_assignment_raw.pdf'), width=w, height=h)
# png(paste0(imgpref,'cancer_assignment.png'), res=450, units='in', width=w, height=h)
# par(mar=c(4,4,0,0))
par(mar=c(4,4,0,0))
image(t(plt.mat), axes=F, useRaster=T, col=hcl.colors(12, "YlOrRd"))
box(lwd=.5)
xat = seq(0,1, length.out=ncol(ca.grp))
yat = seq(0,1, length.out=nrow(ca.grp))
text(y=par()$usr[3] - .005 * diff(par()$usr[3:4]), 
     x=xat, labels=colnames(plt.mat), xpd=TRUE, srt=90, adj=1, cex=.4)
text(x=par()$usr[1] - .005 * diff(par()$usr[1:2]), 
     y=yat, labels=meta[rownames(plt.mat), 'infoline'], xpd=TRUE, adj=1, cex=.25)
# HIGHLIGHT WHERE THE TOP ENRICHMENT SHOULD BE.
gind = sapply(meta[rownames(plt.mat), 'SECONDARY'], function(x){which(colnames(plt.mat) == x)})
points(xat[gind], yat, pch=22, cex=1, lwd=.5)
# Max:
mind = apply(-plt.mat, 1, which.max)
points(xat[mind], yat, pch='*', cex=1, lwd=.25)
dev.off()







# Reorder matrix:
dt <- as.dist(full)
ht <- hclust(dt, method='ward.D')
co <- order.optimal(dt, ht$merge)
cocl = co$order
ht$merge = co$merge
ht$order = cocl
rn <- names(cocl)[cocl]
full = full[rn,rn]
labels = meta[rn, 'GROUP']
faclabels = as.matrix(as.numeric(labels))
colset = as.character(rdcol$COLOR)
# Generate breaks:
NCLUST = 20
breaks = calc.breaks(ht, NCLUST, cocl)
acut <- cutree(ht, NCLUST)[cocl]

# Compute the MDS for observed marks:
mds = list()
c <- cmdscale(as.dist(full))
c <- data.frame(c,id=rownames(c),clust=factor(acut[rownames(c)]))
c = merge(c, meta)
c$COLOR = as.character(c$COLOR)
c$COLOR[c$COLOR == 'white'] <- 'lightgrey'
mds[['Full']] = c
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    mat = obsll[[mark]]
    # Restrict to chosen cells
    cn = rownames(mat)[rownames(mat) %in% cellorder]
    mat = mat[cn, cn]
    c <- cmdscale(as.dist(mat))
    c <- data.frame(c,id=rownames(c),clust=factor(acut[rownames(c)]))
    c = merge(c, meta)
    c$COLOR = as.character(c$COLOR)
    c$COLOR[c$COLOR == 'white'] <- 'lightgrey'
    mds[[mark]] = c
}












