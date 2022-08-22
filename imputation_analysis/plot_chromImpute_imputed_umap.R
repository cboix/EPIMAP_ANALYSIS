#!/usr/bin/R
# ---------------------------------------
# Plot the imputed + observed UMAP
# Plot similar MDS
# ---------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'auxiliary_function_general_repel.R'))
library(ggplot2)
library(viridis)
library(ggrepel)
library(dendextend)
library(uwot)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(ggpubr)
library(igraph)
options(scipen=45) # So we dont get issues writing integers into bedfiles

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    fixedstate = TRUE
    nregions = 0 
    dataset = 'spearman'
    print(paste("No arguments supplied.", 
                "Defaulting to loading nregions =", 
                nregions, 'and dataset', dataset))
} else {        
    fixedstate = as.logical(args[1])
    nregions = as.numeric(args[2])
    dataset = args[3]
}

# Load specific distance matrices:
if (fixedstate){
    commandArgs <- function(trailingOnly=TRUE){ c(nregions, dataset) }
    source(paste0(bindir, 'load_region_distance_matrices.R'))
    setprefix = paste('region',nregions, dataset,'distances_', sep="_")
} else { 
    source(paste0(bindir, 'load_distance_matrices.R'))
    setprefix = 'distances_'
    metric = 'genome-wide correlation'
}

today <- format(Sys.time(), "%m%d%y")
imgdir = paste0(img, "imp_distance/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, setprefix)

parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}

# -------------------------
# Multidimensional scaling:
# -------------------------
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

png(paste0(imgpref,'mds_marks_col.png'),res=300,units='in',width=14,height=5)
layout(matrix(c(1,1,1,1,1,1, 2:16),3,7), heights=c(2,2,2), widths=c(3,3,2,2,2,2,2), TRUE)
par(mar=c(0.25, 0.25, 2, 0.25))
c = mds[['Full']]
c = c[order(c$GROUP, decreasing=TRUE),]
plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='')
# box(lwd=0.5)
# Add a few text (not NONE)
labdf = c[c$GROUP != 'NONE',]
labdf = labdf[order(labdf$X2, labdf$X1),]
sq = seq(0, nrow(labdf), 10)
labdf = labdf[sq,]
text(labdf$X1, labdf$X2, col=labdf$COLOR, 
     labels=labdf$infoline, adj=c(.5,1.3))
mtext('Full', side=3, cex=1)
for (i in 1:(length(marks) -1)){ 
    par(mar=c(0.25, 0.25, 2, 0.25))
    mark = marks[i]
    c = mds[[mark]]
    c = c[order(c$GROUP, decreasing=TRUE),]
    plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='')
    mtext(mark, side=3, cex=1)
    # box(lwd=0.5)
}
dev.off()


# MDS from IMPUTED DATA:
# NOTE: PREVIOUSLY USED DISTANCE ON FULL DATA. CHECK VALIDITY IN MDS FRAMEWORK
imputedmds = list()
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    mat = ll[[mark]]
    # Restrict to chosen cells
    cn = rownames(mat)[rownames(mat) %in% cellorder]
    mat = mat[cn, cn]
    c <- cmdscale(as.dist(mat))
    c <- data.frame(c,id=rownames(c),clust=factor(acut[rownames(c)]))
    c = merge(c, meta)
    c$COLOR = as.character(c$COLOR)
    c$COLOR[c$COLOR == 'white'] <- 'lightgrey'
    imputedmds[[mark]] = c
}

png(paste0(imgpref,'mds_imputedmarks_col.png'),res=300,units='in',width=14,height=5)
layout(matrix(c(1,1,1,1,1,1, 2:16),3,7), heights=c(2,2,2), widths=c(3,3,2,2,2,2,2), TRUE)
par(mar=c(0.25, 0.25, 2, 0.25))
c = mds[['Full']]
c = c[order(c$GROUP, decreasing=TRUE),]
plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='')
labdf = c[c$GROUP != 'NONE',]
labdf = labdf[order(labdf$X2, labdf$X1),]
sq = seq(0, nrow(labdf), 8)
labdf = labdf[sq,]
text(labdf$X1, labdf$X2, col=labdf$COLOR, 
     labels=labdf$infoline, adj=c(.5,1.3))
mtext('Full', side=3, cex=1)
for (i in 1:(length(marks)-1)){ 
    mark = marks[i]
    par(mar=c(0.25, 0.25, 2, 0.25))
    c = imputedmds[[mark]]
    c = c[order(c$GROUP, decreasing=TRUE),]
    plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='')
    mtext(mark, side=3, cex=1)
}
dev.off()

png(paste0(imgpref,'mds_bothmarks_col.png'),res=300,units='in',width=11,height=6.5)
layout(matrix(c(rep(1,4), 2:17),4,4), heights=c(2,2,2,2), widths=c(8,rep(2,3)), TRUE)
par(mar=c(0.25, 0.25, 2.25, 0.25))
c = mds[['Full']]
c = c[order(c$GROUP, decreasing=TRUE),]
plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='')
# Add a few text (not NONE)
labdf = c[c$GROUP != 'NONE',]
labdf = labdf[order(labdf$X2, labdf$X1),]
sq = seq(0, nrow(labdf), 8)
labdf = labdf[sq,]
text(labdf$X1, labdf$X2, col=labdf$COLOR, 
     labels=labdf$infoline, adj=c(.5,1.3))
mtext('Full', side=3, cex=1.5, line=1)
mtext('Imputed', side=3, cex=.75, col='darkgrey')
for (i in 1:NMAIN){ 
    mark = mainmarks[i]
    c = mds[[mark]]
    c = c[order(c$GROUP, decreasing=TRUE),]
    par(mar=c(0.25, 0.25, 2.25, 0.25))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='', cex=.5)
    mtext(mark, side=3, cex=1, line=1)
    mtext('Observed', side=3, cex=.75, col='darkgrey')
    c = imputedmds[[mark]]
    c = c[order(c$GROUP, decreasing=TRUE),]
    par(mar=c(0.25, 0.25, 1, 0.25))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='', cex=.5)
    mtext('Imputed', side=3, cex=.75, col='darkgrey')
}
dev.off()

png(paste0(imgpref,'mds_bothmarks_col_v2.png'),res=300,units='in',width=14,height=8)
layout(matrix(c(1:4, 13, 5:8, 13, 9:12, 13),3,5, byrow=T) , heights=c(2,2,2), widths=c(rep(2,4),4), TRUE)
for (i in 1:NMAIN){ 
    mark = mainmarks[i]
    c = mds[[mark]]
    c = c[order(c$GROUP, decreasing=TRUE),]
    par(mar=c(0.25, 2.25, 1, 0.25))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='', cex=.75)
    mtext(mark, side=2, cex=1, line=0)
    mtext('Observed', side=3, cex=.75, col='darkgrey')
    c = imputedmds[[mark]]
    c = c[order(c$GROUP, decreasing=TRUE),]
    par(mar=c(0.25, 0.25, 1, 0.25))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='', cex=.75)
    mtext('Imputed', side=3, cex=.75, col='darkgrey')
}
par(mar=c(0.25, 0.25, 2.25, 0.25))
c = mds[['Full']]
c = c[order(c$GROUP, decreasing=TRUE),]
plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='')
# Add a few text (not NONE)
labdf = c[c$GROUP != 'NONE',]
labdf = labdf[order(labdf$X2, labdf$X1),]
sq = seq(0, nrow(labdf), 8)
labdf = labdf[sq,]
text(labdf$X1, labdf$X2, col=labdf$COLOR, 
     labels=labdf$infoline, cex=.8, adj=c(.5,1.3))
mtext('Full', side=3, cex=1.5, line=1)
mtext('Imputed', side=3, cex=.75, col='darkgrey')
dev.off()


# -------------------------
# UMAP embeddings:
# -------------------------
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

run.umap = function(mat, acut, meta, nn=150, mdist=0.25, learning_rate=0.5, separate.io=FALSE){
    diag(mat) = 0
    u = umap(as.dist(mat), n_neighbors=nn, min_dist=mdist, verbose=F,
             repulsion_strength=.25)
    if (separate.io){
        mdf = t(sapply(rownames(mat), function(x){strsplit(x, "_")[[1]]}))
        udf <- data.frame(X1=u[,1],X2=u[,2], dataset=mdf[,2],
                          id=mdf[,1],clust=factor(acut[mdf[,1]]))
    } else {
        udf <- data.frame(X1=u[,1],X2=u[,2],
                          id=rownames(mat),clust=factor(acut[rownames(mat)]))
    }
    udf = merge(udf, meta)
    udf$COLOR = as.character(udf$COLOR)
    udf$COLOR[udf$COLOR == 'white'] <- 'lightgrey'
    return(udf)
}

umap = list()
umap[['Full']] = run.umap(full, acut, meta)
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    mat = obsll[[mark]]
    # Restrict to chosen cells
    cn = rownames(mat)[rownames(mat) %in% cellorder]
    mat = mat[cn, cn]
    umap[[mark]] = run.umap(mat, acut, meta, nn=25)
}

# UMAP from Imputed data:
imputedumap = list()
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    mat = ll[[mark]]
    cn = rownames(mat)[rownames(mat) %in% cellorder]
    mat = mat[cn, cn]
    imputedumap[[mark]] = run.umap(mat, acut, meta, nn=200)
}

png(paste0(imgpref,'umap_bothmarks_col_v2.png'),res=300,units='in',width=14,height=8)
layout(matrix(c(1:4, 13, 5:8, 13, 9:12, 13),3,5, byrow=T) , heights=c(2,2,2), widths=c(rep(2,4),4), TRUE)
for (i in 1:NMAIN){ 
    mark = mainmarks[i]
    c = umap[[mark]]
    c = c[order(c$GROUP, decreasing=TRUE),]
    par(mar=c(0.25, 2.25, 1, 0.25))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='', cex=.75)
    mtext(mark, side=2, cex=1, line=0)
    mtext('Observed', side=3, cex=.75, col='darkgrey')
    c = imputedumap[[mark]]
    c = c[order(c$GROUP, decreasing=TRUE),]
    par(mar=c(0.25, 0.25, 1, 0.25))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='', cex=.75)
    mtext('Imputed', side=3, cex=.75, col='darkgrey')
}
par(mar=c(0.25, 0.25, 2.25, 0.25))
c = umap[['Full']]
c = c[order(c$GROUP, decreasing=TRUE),]
plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='')
# Add a few text (not NONE)
labdf = c[c$GROUP != 'NONE',]
labdf = labdf[order(labdf$X2, labdf$X1),]
sq = seq(0, nrow(labdf), 8)
labdf = labdf[sq,]
text(labdf$X1, labdf$X2, col=labdf$COLOR, 
     labels=labdf$infoline, cex=.8, adj=c(.5,1.3))
mtext('Full', side=3, cex=1.5, line=1)
mtext('Imputed', side=3, cex=.75, col='darkgrey')
dev.off()


# ---------------------------------
# Run for the same samples in both:
# For reviews
# ---------------------------------
boumap = list()
biumap = list()
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    print(mark)
    omat = obsll[[mark]]
    imat = ll[[mark]]
    # Restrict to chosen cells
    cn = rownames(omat)[rownames(omat) %in% cellorder]
    omat = omat[cn, cn]
    icn = rownames(imat)[rownames(imat) %in% cellorder]
    imat = imat[icn, icn]
    bcn = icn[icn %in% cn]
    names(bcn) = NULL
    bomat = omat[bcn, bcn]
    bimat = imat[bcn, bcn]
    boumap[[mark]] = run.umap(bomat, acut, meta, nn=25)
    biumap[[mark]] = run.umap(bimat, acut, meta, nn=25)
}

png(paste0(imgpref,'umap_samesamp_bothmarks_col.png'),res=300,units='in',width=8,height=14)
# layout(matrix(c(1:4, 13, 5:8, 13, 9:12, 13),3,5, byrow=T) , heights=c(2,2,2), widths=c(rep(2,4),4), TRUE)
layout(matrix(c(1:28),7,4, byrow=T) , heights=rep(2,7), widths=rep(2,4), TRUE)
sp=0
# for (i in 1:NMAIN){ 
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    # Observed:
    # mark = mainmarks[i]
    c = boumap[[mark]]
    c = c[order(c$GROUP, decreasing=TRUE),]
    par(mar=c(sp, 2.25, 1, sp))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='', cex=.75)
    mtext(mark, side=2, cex=1, line=0)
    mtext('Observed', side=3, cex=.75, col='darkgrey')
    box()
    # Imputed:
    c = biumap[[mark]]
    c = c[order(c$GROUP, decreasing=TRUE),]
    par(mar=c(sp, sp, 1, sp))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='', cex=.75)
    mtext('Imputed', side=3, cex=.75, col='darkgrey')
    box()
}
dev.off()


# ------------------------------------------------------------
# Calculate simple clustering-coefficient/homogeneity metrics:
# ------------------------------------------------------------
# 1. TODO: Clustering coefficient (using the N nearest neighbors)
# 2. Within-group distance vs. Out-group distance?
# Find samples with similar properties:
llist = sapply(cellorder, function(x){
                   cind = which((meta[cellorder,'infoline'] == meta[x,'infoline']) * 
                                (meta[cellorder,'lifestage'] == meta[x,'lifestage']) *
                                (meta[cellorder,'GROUP'] == meta[x,'GROUP']) == 1)
                   cellorder[cind] } )
llist = llist[lapply(llist, length) > 2]
llist = unique(llist)

rlist = lapply(llist, function(x){ sample(cellorder, length(x))})
meta$GROUP = factor(meta$GROUP, levels=odf$GROUP)

# Plot ordered by group:
png(paste0(imgpref, 'matrices_grpord_allmarks.png'),res=300,units='in',width=13,height=2.2)
layout(matrix(1:(13 * 3), nrow=3, byrow=F), heights=c(.2, 1,1), widths=rep(1,13), TRUE)
for (i in 1:(length(marks) -1)){
    # Load both observed and imputed datasets:
    mark = marks[i]
    print(mark)
    omat = obsll[[mark]]
    imat = ll[[mark]]
    # Restrict to chosen cells:
    cn = rownames(omat)[rownames(omat) %in% cellorder]
    omat = omat[cn, cn]
    icn = rownames(imat)[rownames(imat) %in% cellorder]
    imat = imat[icn, icn]
    # Also calculate the same exact restricted matrix:
    bcn = icn[icn %in% cn]
    names(bcn) = NULL
    bomat = omat[bcn, bcn]
    bimat = imat[bcn, bcn]
    bicn = rownames(reord(bimat))
    # Plot, ordered by group:
    grpord = bicn[order(meta[bicn,'GROUP'])]
    # Plot 2x2 ordering with the colors on the side...
    sp=0.1
    par(mar=c(sp, sp,1,sp))
    meta.image(metamat[grpord,5,drop=F], colvals=colvals, cex=0, horiz=F, useRaster=TRUE)
    mtext(paste(mark), cex=.75)
    par(mar=c(sp, sp,sp,sp))
    image(bomat[grpord, grpord], useRaster=T, axes=F, col=colspec)
    image(bimat[grpord, grpord], useRaster=T, axes=F, col=colspec)
}
dev.off()



ccdf = c()
all.wdf = c()
all.samedf = c()
all.randdf = c()
all.wodf = c()
for (i in 1:(length(marks) -1)){
    # Load both observed and imputed datasets:
    mark = marks[i]
    print(mark)
    omat = obsll[[mark]]
    imat = ll[[mark]]
    # Restrict to chosen cells:
    cn = rownames(omat)[rownames(omat) %in% cellorder]
    omat = omat[cn, cn]
    icn = rownames(imat)[rownames(imat) %in% cellorder]
    imat = imat[icn, icn]
    # Also calculate the same exact restricted matrix:
    bcn = icn[icn %in% cn]
    names(bcn) = NULL
    bomat = omat[bcn, bcn]
    bimat = imat[bcn, bcn]
    bocn = rownames(reord(bomat))
    bicn = rownames(reord(bimat))

    # Plot 2x2 ordering with the colors on the side...
    png(paste0(imgpref, 'matrices_samerows_impobs_', mark,'.png'),res=300,units='in',width=5.5,height=5)
    layout(matrix(1:6, ncol=3, byrow=T), heights=c(1,1), widths=c(.15,1,1), TRUE)
    sp=0.1
    par(mar=c(sp, sp,1,sp))
    meta.image(metamat[bocn,5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
    image(bomat[bocn, bocn], useRaster=T, axes=F, col=colspec)
    mtext(paste("Observed", mark, "(ordered by Observed)"), cex=.75)
    image(bimat[bocn, bocn], useRaster=T, axes=F, col=colspec)
    mtext(paste("Imputed", mark, "(ordered by Observed)"), cex=.75)
    # Ord by imputed
    meta.image(metamat[bicn,5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
    image(bomat[bicn, bicn], useRaster=T, axes=F, col=colspec)
    mtext(paste("Observed", mark, "(ordered by Imputed)"), cex=.75)
    image(bimat[bicn, bicn], useRaster=T, axes=F, col=colspec)
    mtext(paste("Imputed", mark, "(ordered by Imputed)"), cex=.75)
    dev.off()

    # Plot, ordered by group:
    grpord = bicn[order(meta[bicn,'GROUP'])]
    # Plot 2x2 ordering with the colors on the side...
    png(paste0(imgpref, 'matrices_grpord_', mark,'.png'),res=300,units='in',width=5.5,height=2.5)
    layout(matrix(1:3, ncol=3, byrow=T), heights=c(1), widths=c(.15,1,1), TRUE)
    sp=0.1
    par(mar=c(sp, sp,1,sp))
    meta.image(metamat[grpord,5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
    image(bomat[grpord, grpord], useRaster=T, axes=F, col=colspec)
    mtext(paste("Observed", mark), cex=.75)
    image(bimat[grpord, grpord], useRaster=T, axes=F, col=colspec)
    mtext(paste("Imputed", mark), cex=.75)
    dev.off()


    # Calculate homogeneity statistics:
    # On same set, calc avg. distance between members of the same set
    tform = make.tform(metamat[bcn,'group'], norm=TRUE, u=odf$GROUP)
    # TODO: Mean distance and variance in distances?
    bomat.norm = bomat - mean(bomat)
    bimat.norm = bimat - mean(bimat)
    # gomat = t(tform) %*% (bomat.norm %*% tform)
    # gimat = t(tform) %*% (bimat.norm %*% tform)
    gomat = t(tform) %*% (bomat %*% tform)
    gimat = t(tform) %*% (bimat %*% tform)

    # Only within dist --> (TODO) show fit:
    wdf = data.frame(obs=diag(gomat), imp=diag(gimat), mark=mark)
    wdf$GROUP = rownames(wdf)
    wdf =merge(wdf, odf)
    wdf = rbind(all.wdf, wdf)
    fit = lm(imp ~ obs, wdf)
    print(summary(fit))

    png(paste0(imgpref, 'avg_group_dist_impobs_', mark,'.png'),res=300,units='in',width=8,height=5)
    layout(matrix(1:2, ncol=2))
    par(mar=c(4,4,2,1))
    plot(wdf$obs,wdf$imp, pch=19, xlab='Avg. within dist (Observed)', ylab='Avg. within dist (Imputed)',
         col=odf$COLOR, xlim=c(0,1), ylim=c(0,1), bty='n')
    mtext(paste(mark, '(Only within-group distances)'))
    # Multiple comparisons of between/within dist:
    par(mar=c(4,4,2,1))
    plot(gomat, gimat, pch=19, xlab='Avg. dist (Observed)', ylab='Avg. dist (Imputed)', 
         xlim=c(0,1), ylim=c(0,1), bty='n')
    abline(0,1)
    mtext(paste(mark, '(All cross-group distances)'))
    dev.off()

    # All comparisons:
    bodf = data.frame(gomat)
    bodf$GROUP = rownames(bodf)
    bodf = gather(bodf, group2, dist, -GROUP)
    bodf = merge(bodf, odf)
    bodf$type = 'Observed'
    # 
    bidf = data.frame(gimat)
    bidf$GROUP = rownames(bidf)
    bidf = gather(bidf, group2, dist, -GROUP)
    bidf = merge(bidf, odf)
    bidf$type = 'Imputed'
    bdf = rbind(bodf, bidf)
    bdf = bdf[!is.na(bdf$dist),]

    gplot = ggplot(bdf, aes(group2, dist, color=group2)) + 
        geom_point() + 
        facet_grid(type ~ GROUP, scales='free_y') + 
        scale_color_manual(values=colvals$group) + 
        labs('Average cross-cluster distance') +
        theme_pubclean()  + 
        theme(axis.text.x  = element_blank(),
              axis.ticks.x  = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              # axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA, size=1))
    ggsave(paste0(imgpref, 'avgdist_crossgroup', mark, '.png'), gplot, dpi=450, units='in', width=15, height=5)


    # Calc some within-group vs out-group metrics:
    wometrics = sapply(odf$GROUP, function(x){
                           samp = bcn[metamat[bcn, 'group'] == x]
                           nosamp = bcn[metamat[bcn, 'group'] != x]
                           # Eval in-group/out-group in both types:
                           somat = bomat[samp, samp]
                           simat = bimat[samp, samp]
                           with.obs = mean(somat[lower.tri(somat)])
                           with.imp = mean(simat[lower.tri(simat)])
                           wv.obs = var(somat[lower.tri(somat)])
                           wv.imp = var(simat[lower.tri(simat)])
                           nomat = bomat[samp, nosamp]
                           nimat = bimat[samp, nosamp]
                           out.obs = mean(nomat)
                           out.imp = mean(nimat)
                           ov.obs = var(c(nomat))
                           ov.imp = var(c(nimat))
                           return(c(with.obs / out.obs, with.imp / out.imp,
                                    wv.obs / ov.obs, wv.imp / ov.imp))
              })


    png(paste0(imgpref, 'with_out_mean_impobs_', mark,'.png'),res=300,units='in',width=5.5,height=5)
    par(mar=c(4,4,1.5,1))
    rn = max(wometrics[1:2,], na.rm=T) * 1.1
    plot(wometrics[1,], wometrics[2,], pch=19, col=odf$COLOR, 
         xlab='Avg. dist in-group / out-group (Observed)', 
         ylab='Avg. dist in-group / out-group (Imputed)', 
         xlim=c(0,rn), ylim=c(0,rn), bty='n')
    abline(0,1)
    mtext(paste(mark,'(in/out-group mean distance)'))
    dev.off()

    png(paste0(imgpref, 'with_out_var_impobs_', mark,'.png'),res=300,units='in',width=5.5,height=5)
    rn = max(wometrics[3:4,], na.rm=T) * 1.1
    par(mar=c(4,4,1.5,1))
    plot(wometrics[3,], wometrics[4,], pch=19, col=odf$COLOR, 
         xlab='Avg. var in-group / out-group (Observed)', 
         ylab='Avg. var in-group / out-group (Imputed)', 
         xlim=c(0,rn), ylim=c(0,rn), bty='n')
    abline(0,1)
    axis(1)
    axis(2)
    mtext(paste(mark,'(in/out-group variance in distance)'))
    dev.off()

    wodf = data.frame(t(wometrics), mark=mark)
    names(wodf)[1:4] = c('wm','om','wv','ov')
    wodf$GROUP = rownames(wodf)
    wodf = merge(wodf, odf)
    all.wodf = rbind(all.wodf, wodf)

    # Calculate homogeneity statistics:

    # Clustering coefficient:
    NTOP = 10
    # TODO: Normalize matrices?:
    # bonorm = sweep(bomat,1,apply(bomat,1,mean), '-')
    # bonorm = sweep(bonorm,2,apply(bomat,1,mean), '-')
    # binorm = sweep(bimat,1,apply(bimat,1,mean), '-')
    # binorm = sweep(binorm,2,apply(bimat,1,mean), '-')

    # Make adjacency:
    # oadj = apply(bonorm, 1, function(x){
    oadj = apply(bomat, 1, function(x){
                     thresh = sort(x)[NTOP]
                     1 * (x <= thresh) })

    # iadj = apply(binorm, 1, function(x){
    iadj = apply(bimat, 1, function(x){
                     thresh = sort(x)[NTOP]
                     1 * (x <= thresh) })

    og = graph_from_adjacency_matrix(oadj)
    ig = graph_from_adjacency_matrix(iadj)

    occ = transitivity(og)
    icc = transitivity(ig)
    ccdf = rbind(ccdf, data.frame(obs=occ, imp=icc, mark=mark))

    # locc = transitivity(og, type='local')
    # licc = transitivity(ig, type='local')
    # plot(locc, licc, pch=19, col=meta[vertex_attr(og)$name,'COLOR'])

    # Look at similar samples:
    slist = sapply(1:length(llist), function(i){
                       x = llist[[i]]
                       samp = bcn[bcn %in% x]
                       nosamp = bcn[!(bcn %in% samp)]
                       if (length(samp) > 1){
                           somat = bomat[samp, samp]
                           simat = bimat[samp, samp]
                           with.obs = mean(somat[lower.tri(somat)])
                           with.imp = mean(simat[lower.tri(simat)])
                           nomat = bomat[samp, nosamp]
                           nimat = bimat[samp, nosamp]
                           out.obs = mean(nomat)
                           out.imp = mean(nimat)
                           return(c(i,meta[x[1], c('infoline')],
                                    with.obs / out.obs, with.imp / out.imp)) }})
    samedf = data.frame(do.call(rbind, slist))
    samedf$mark = mark
    all.samedf = rbind(all.samedf, samedf)

    rslist = sapply(1:length(rlist), function(i){
                       x = rlist[[i]]
                       samp = bcn[bcn %in% x]
                       nosamp = bcn[!(bcn %in% samp)]
                       if (length(samp) > 1){
                           somat = bomat[samp, samp]
                           simat = bimat[samp, samp]
                           with.obs = mean(somat[lower.tri(somat)])
                           with.imp = mean(simat[lower.tri(simat)])
                           nomat = bomat[samp, nosamp]
                           nimat = bimat[samp, nosamp]
                           out.obs = mean(nomat)
                           out.imp = mean(nimat)
                           return(c(i,meta[x[1], c('infoline')],
                                    with.obs / out.obs, with.imp / out.imp)) }})
    randdf = data.frame(do.call(rbind, rslist))
    randdf$mark = mark
    all.randdf = rbind(all.randdf, randdf)
}



# Clustering coefficient plots:
gplot = ggplot(ccdf, aes(obs, imp, label=mark)) + 
    geom_point() + 
    geom_text_repel() + 
    xlim(0,1) + ylim(0,1) + 
    geom_abline(intercept=0,slope=1) + 
    labs(x='Observed', y='Imputed', title='Clustering Coefficient (graph with N=25)') + 
    theme_pubr()
ggsave(paste0(imgpref, 'transitivity_impobs_allmarks.png'), gplot,dpi=300,units='in',width=5,height=5)

ccdf = ccdf[order(ccdf$imp/ccdf$obs),]
ccdf$mark = factor(ccdf$mark, levels=ccdf$mark)

gplot = ggplot(ccdf, aes(mark, imp/obs)) + 
    geom_bar(stat='identity', fill='grey75') + 
    geom_hline(yintercept=1) + 
    labs(x='Mark/Assay', y='Imputed/Observed Coefficient', title='Ratio of Clust. Coeff. (graph with N=25)') + 
    theme_pubr() + theme(axis.text.x=element_text(angle=90, hjust=1)) 
ggsave(paste0(imgpref, 'transitivity_ratio_impobs_allmarks.png'), gplot,dpi=300,units='in',width=5,height=5)


# Plot the within mark/without etc:
gplot = ggplot(all.wodf, aes(wm, om, color=GROUP)) + 
    geom_point() + 
    geom_abline() + 
    scale_color_manual(values=colvals$group) +
    labs(x='Avg. dist in-group / out-group (Observed)', y='Avg. dist in-group / out-group (Imputed)') +
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'inoutgroup_mean.png'), gplot,dpi=300,units='in',width=5,height=5)


# Plot the within mark/without etc:
gplot = ggplot(all.wodf, aes(wv, ov, color=GROUP)) + 
    geom_point() + 
    geom_abline() + 
    scale_color_manual(values=colvals$group) +
    labs(x='Var. dist. in-group / out-group (Observed)', y='Var. dist. in-group / out-group (Imputed)') +
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'inoutgroup_var.png'), gplot,dpi=300,units='in',width=5,height=5)


names(all.samedf) = c('num','infoline','obs','imp','mark')
all.samedf = merge(all.samedf, unique(meta[,c('GROUP','infoline')])) 
all.samedf = merge(all.samedf, odf) 
all.samedf$obs = as.numeric(as.character(all.samedf$obs))
all.samedf$imp = as.numeric(as.character(all.samedf$imp))

names(all.randdf) = c('num','infoline','obs','imp','mark')
all.randdf = merge(all.randdf, unique(meta[,c('GROUP','infoline')])) 
all.randdf = merge(all.randdf, odf) 
all.randdf$obs = as.numeric(as.character(all.randdf$obs))
all.randdf$imp = as.numeric(as.character(all.randdf$imp))


# Plot close biological samples (TODO: compare to rand samples or samples within same group)
mx = max(c(all.samedf$obs, all.samedf$imp))
gplot = ggplot(all.samedf, aes(obs, imp, color=GROUP)) + 
    geom_point() + 
    geom_abline() + 
    xlim(0, mx) + ylim(0,mx)  +
    scale_color_manual(values=colvals$group) +
    labs(x='Avg. dist in/out-group (Observed)', y='Avg. dist in/out-group (Imputed)') +
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'inout_sametype_mean.png'), gplot,dpi=300,units='in',width=5,height=5)

# Comparison if we do random samples of same size:
mx = max(c(all.randdf$obs, all.randdf$imp))
gplot = ggplot(all.randdf, aes(obs, imp)) + 
    geom_point() + 
    geom_abline() + 
    xlim(0, mx) + ylim(0,mx)  +
    scale_color_manual(values=colvals$group) +
    labs(x='Avg. dist in/out-group (Observed)', y='Avg. dist in/out-group (Imputed)') +
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'inout_randtype_mean.png'), gplot,dpi=300,units='in',width=5,height=5)


# Run for the same samples in both:
matchumap = list()
NNVAL = 50
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    print(mark)
    omat = obsll[[mark]]
    imat = ll[[mark]]
    # Restrict to chosen cells
    cn = rownames(omat)[rownames(omat) %in% cellorder]
    icn = rownames(imat)[rownames(imat) %in% cellorder]
    bcn = icn[icn %in% cn]
    names(bcn) = NULL
    mat = full.ll[[mark]]
    samp = t(sapply(rownames(mat),
                    function(x){strsplit(x, "_")[[1]][1]}))
    cn = rownames(mat)[samp %in% bcn]
    mat = mat[cn, cn]
    set.seed(1)
    matchumap[[mark]] = run.umap(mat, acut, meta, nn=NNVAL, separate.io=TRUE)
}


# Plot the observed/imputed + links:
png(paste0(imgpref,'umap_links_allmarks_n',NNVAL,'.png'),res=300,units='in',width=9,height=9)
layout(matrix(c(1:16),4,4, byrow=T) , heights=rep(2,4), widths=rep(2,4), TRUE)
sp=0
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    udf = matchumap[[mark]]
    udf$COLOR = as.character(udf$COLOR)
    udf$COLOR[udf$COLOR == 'white'] <- 'lightgrey'
    udf$pch = 19
    udf$pch[udf$dataset == 'imp'] = 1
    # Separate out matched segments:
    count = aggregate(dataset ~ id, udf, function(x){length(unique(x))})
    keepct = as.character(count$id[count$dataset == 2])
    idf = merge(data.frame(id=keepct, dataset="imp"), udf[,c('id','dataset','X1','X2','COLOR')])
    colnames(idf)[3:4] = c('I1','I2')
    obdf = merge(data.frame(id=keepct, dataset="obs"), udf[,c('id','dataset','X1','X2','COLOR')])
    segdf = merge(obdf[,c('id','X1','X2','COLOR')], idf[,c('id','I1','I2')])
    c = udf
    c = c[order(c$GROUP, decreasing=TRUE),]
    # Print out ones with most len:
    segdf$len = sqrt((segdf$X1 - segdf$I1)^2 + (segdf$X2 - segdf$I2)^2)
    segdf = merge(segdf, meta)
    ssdf = segdf[segdf$len > 1,]
    xlim = c(min(c$X1), max(c$X1))
    ylim = c(min(c$X2), max(c$X2))
    # Make actual plot
    par(mar=c(sp,sp,1.5,sp))
    plot(c$X1, c$X2, col=c$COLOR, pch=c$pch, ylab='UMAP-2', xlab='UMAP-1',
         cex=.75, xlim=xlim, ylim=ylim, axes=F)
    segments(segdf$X1,segdf$X2, segdf$I1, segdf$I2, col=segdf$COLOR, lwd=.5)
    box(lwd=.5)
    # mtext(paste0('UMAP on ', metric,' of ', mark, '\n Imputed (open) and Observed (closed)'), side=3, line=1)
    mtext(paste0(mark))
    if (nrow(ssdf) > 0){
        # text(ssdf$X1, ssdf$X2, labels=paste0(ssdf$id, "(", ssdf$GROUP,")\n", ssdf$infoline), col=ssdf$COLOR, cex=.5)
    }
}
plot(1,1, type='n', axes=F)
legend('center', c('Observed','Imputed'), pch=c(19,1), 
       col='black', pt.cex=2, cex=1.5, bty='n')
dev.off()

# Plot the colors, no links:
png(paste0(imgpref,'umap_links_allmarks_group_n',NNVAL,'.png'),res=300,units='in',width=7,height=2.2)
layout(matrix(c(1:14),2,7, byrow=T) , heights=rep(1.1,2), widths=rep(1,7), TRUE)
sp=0.05
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    udf = matchumap[[mark]]
    udf$COLOR = as.character(udf$COLOR)
    udf$COLOR[udf$COLOR == 'white'] <- 'lightgrey'
    udf$pch = 19
    udf$pch[udf$dataset == 'imp'] = 1
    # Separate out matched segments:
    count = aggregate(dataset ~ id, udf, function(x){length(unique(x))})
    keepct = as.character(count$id[count$dataset == 2])
    idf = merge(data.frame(id=keepct, dataset="imp"), udf[,c('id','dataset','X1','X2','COLOR')])
    colnames(idf)[3:4] = c('I1','I2')
    obdf = merge(data.frame(id=keepct, dataset="obs"), udf[,c('id','dataset','X1','X2','COLOR')])
    segdf = merge(obdf[,c('id','X1','X2','COLOR')], idf[,c('id','I1','I2')])
    c = udf
    ind = 1:nrow(c)
    ind = sample(ind, replace=FALSE)
    c = c[ind,]
    # Print out ones with most len:
    segdf$len = sqrt((segdf$X1 - segdf$I1)^2 + (segdf$X2 - segdf$I2)^2)
    segdf = merge(segdf, meta)
    ssdf = segdf[segdf$len > 1,]
    xlim = c(min(c$X1), max(c$X1))
    ylim = c(min(c$X2), max(c$X2))
    # Make actual plot
    par(mar=c(sp,sp,1.5,sp))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, ylab='UMAP-2', xlab='UMAP-1',
         cex=.4, xlim=xlim, ylim=ylim, axes=F)
    rect(xleft=parpos(1,0), xright=parpos(1,-1),
         ybottom=parpos(2,-1), ytop=parpos(2,-1.3), col='grey80', xpd=TRUE, border=NA)
    mtext(paste0(mark))
}
plot(1,1, type='n', axes=F)
legend('top', odf$GROUP, pch=c(19,19), inset=c(0,-.25), ncol=2, xpd=TRUE,
       col=odf$COLOR, pt.cex=.75, cex=.45, bty='n')
dev.off()


# Plot the observed/imputed + links:
png(paste0(imgpref,'umap_links_allmarks_impobs_n',NNVAL,'.png'),res=300,units='in',width=7,height=2.2)
layout(matrix(c(1:14),2,7, byrow=T) , heights=rep(1.1,2), widths=rep(1,7), TRUE)
sp=0.05
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    udf = matchumap[[mark]]
    udf$COLOR = as.character(udf$COLOR)
    udf$COLOR[udf$COLOR == 'white'] <- 'lightgrey'
    udf$pch = 19
    udf$pch[udf$dataset == 'imp'] = 1
    # Separate out matched segments:
    count = aggregate(dataset ~ id, udf, function(x){length(unique(x))})
    keepct = as.character(count$id[count$dataset == 2])
    idf = merge(data.frame(id=keepct, dataset="imp"), udf[,c('id','dataset','X1','X2','COLOR')])
    colnames(idf)[3:4] = c('I1','I2')
    obdf = merge(data.frame(id=keepct, dataset="obs"), udf[,c('id','dataset','X1','X2','COLOR')])
    segdf = merge(obdf[,c('id','X1','X2','COLOR')], idf[,c('id','I1','I2')])
    c = udf
    ind = 1:nrow(c)
    ind = sample(ind, replace=FALSE)
    c = c[ind,]
    # Print out ones with most len:
    segdf$len = sqrt((segdf$X1 - segdf$I1)^2 + (segdf$X2 - segdf$I2)^2)
    segdf = merge(segdf, meta)
    ssdf = segdf[segdf$len > 1,]
    xlim = c(min(c$X1), max(c$X1))
    ylim = c(min(c$X2), max(c$X2))
    # Make actual plot
    c$COLOR = ifelse(c$dataset == 'imp','indianred','royalblue')
    par(mar=c(sp,sp,1.5,sp))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, ylab='UMAP-2', xlab='UMAP-1',
         cex=.4, xlim=xlim, ylim=ylim, axes=F)
    rect(xleft=parpos(1,0), xright=parpos(1,-1),
         ybottom=parpos(2,-1), ytop=parpos(2,-1.3), col='grey80', xpd=TRUE, border=NA)
    mtext(paste0(mark))
}
plot(1,1, type='n', axes=F)
legend('top', c('Observed','Imputed'), pch=c(19,19), 
       col=c('royalblue','indianred'), pt.cex=1.5, cex=1, bty='n')
dev.off()



# --------------------
# UMAP from full data:
# --------------------
NNVAL=250
fullumap = list()
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    print(mark)
    mat = full.ll[[mark]]
    # Restrict to chosen cells:
    samp = t(sapply(rownames(mat),
                    function(x){strsplit(x, "_")[[1]][1]}))
    cn = rownames(mat)[samp %in% cellorder]
    mat = mat[cn, cn]
    set.seed(1)
    fullumap[[mark]] = run.umap(mat, acut, meta, nn=NNVAL, separate.io=TRUE)
}

png(paste0(imgpref,'umap_samerun_bothmarks_col_n',NNVAL,'.png'),res=300,units='in',width=10,height=7)
layout(matrix(c(1:12),3,4, byrow=T) , heights=c(2,2,2), widths=c(rep(2,4)), TRUE)
for (i in 1:NMAIN){ 
    mark = mainmarks[i]
    df = fullumap[[mark]]
    df = df[order(df$GROUP, decreasing=TRUE),]
    c = df[df$dataset == 'obs',]
    # Whether to put header or not:
    if (i %% 6 %in% c(1,2)){ lh = TRUE } else { lh = FALSE }
    # Observed:
    par(mar=c(0.25, 2.25, 0.25 + .8 * lh, 0.0))
    # plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='', cex=.75)
    plot(c$X1, c$X2, col=c$COLOR, pch=19, yaxt='n', xaxt='n', ylab='', xlab='', cex=.5)
    mtext(mark, side=2, cex=1.4, line=0)
    if(lh) {mtext('Observed', side=3, cex=1, col='black')}
    # Imputed
    c = df[df$dataset == 'imp',]
    par(mar=c(0.25, 0.0, 0.25 + .8 * lh, 0.25))
    # plot(c$X1, c$X2, col=c$COLOR, pch=19, axes=F, ylab='', xlab='', cex=.75)
    plot(c$X1, c$X2, col=c$COLOR, pch=19, yaxt='n', xaxt='n', ylab='', xlab='', cex=.5)
    if (lh){ mtext('Imputed', side=3, cex=1, col='black') }
}
dev.off()

png(paste0(imgpref,'umap_samerun_bothmarks_col_t12_n',NNVAL,'.png'),res=300,units='in',width=15,height=5)
lmat = cbind(matrix(1:6, 3, 2, byrow=T),
             matrix(1:6, 3, 2, byrow=T) + 6,
             matrix(1:6, 3, 2, byrow=T) + 12,
             matrix(1:6, 3, 2, byrow=T) + 18,
             matrix(1:6, 3, 2, byrow=T) + 24)
layout(lmat, heights=c(2.25,2,2), widths=c(rep(2,10)), TRUE)
for (i in 1:NT12){ 
    mark = t12marks[i]
    df = fullumap[[mark]]
    df = df[order(df$GROUP, decreasing=TRUE),]
    c = df[df$dataset == 'obs',]
    # Whether to put header or not:
    if (i %% 3 == 1){ lh = TRUE } else { lh = FALSE }
    # Observed:
    par(mar=c(0.25, 2.25, 0.25 + 1 * lh, 0.0))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, yaxt='n', xaxt='n', ylab='', xlab='', cex=.5)
    mtext(mark, side=2, cex=1.4, line=0)
    if(lh) {mtext('Observed', side=3, cex=1, col='black')}
    # Imputed
    c = df[df$dataset == 'imp',]
    par(mar=c(0.25, 0.0, 0.25 + 1 * lh, 0.25))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, yaxt='n', xaxt='n', ylab='', xlab='', cex=.5)
    if (lh){ mtext('Imputed', side=3, cex=1, col='black') }
}
dev.off()


# TODO MAKE:
# Observed
# Imputed
# New imputed only
# Paired, with links
png(paste0(imgpref,'umap_samerun_bothmarks_col_t12_3panel_n',NNVAL,'.png'),res=300,units='in',width=17,height=6)
lmat = cbind(matrix(1:12, 4, 3, byrow=T),
             matrix(1:12, 4, 3, byrow=T) + 12,
             matrix(1:12, 4, 3, byrow=T) + 24,
             matrix(1:12, 4, 3, byrow=T) + 36)
layout(lmat, heights=c(2.25,2,2,2), widths=c(rep(2,12)), TRUE)
for (i in 1:NT12){ 
    mark = t12marks[i]
    df = fullumap[[mark]]
    df = df[order(df$GROUP, decreasing=TRUE),]
    xlim = range(df$X1)
    ylim = range(df$X2)
    c = df[df$dataset == 'obs',]
    obsid = as.character(c$id)
    # Whether to put header or not:
    if (i %% 4 == 1){ lh = TRUE } else { lh = FALSE }
    # Observed:
    par(mar=c(0.25, 2.25, 0.25 + 1 * lh, 0.0))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, yaxt='n', xaxt='n', ylab='', xlab='', cex=.5, xlim=xlim, ylim=ylim)
    # plot(c$X1, c$X2, col=c$COLOR, pch=c$lifestage, yaxt='n', xaxt='n', ylab='', xlab='', cex=.5, xlim=xlim, ylim=ylim)
    mtext(mark, side=2, cex=1.4, line=0)
    if(lh) {mtext('Observed', side=3, cex=1, col='black')}
    # Imputed
    c = df[df$dataset == 'imp',]
    par(mar=c(0.25, 0.0, 0.25 + 1 * lh, 0.0))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, yaxt='n', xaxt='n', ylab='', xlab='', cex=.5, xlim=xlim, ylim=ylim)
    if (lh){ mtext('Imputed (all)', side=3, cex=1, col='black') }
    # Only novel
    c = c[!(c$id %in% obsid),]
    par(mar=c(0.25, 0.0, 0.25 + 1 * lh, 0.25))
    plot(c$X1, c$X2, col=c$COLOR, pch=19, yaxt='n', xaxt='n', ylab='', xlab='', cex=.5, xlim=xlim, ylim=ylim)
    if (lh){ mtext('Imputed (novel)', side=3, cex=1, col='black') }
}
dev.off()


plot.largepanel = function(pchval=NULL, colval=NULL, axonly=FALSE, ptonly=FALSE){
    lmat = cbind(matrix(c(rep(1,6),2,3),4,2, byrow=T),
                 3 + cbind(matrix(1:12, 4, 3, byrow=T),
                           matrix(1:12, 4, 3, byrow=T) + 12,
                           matrix(1:12, 4, 3, byrow=T) + 24))
    layout(lmat, heights=c(2.25,2,2,2), widths=c(3,3,rep(2,9)), TRUE)
    i = 1
    if(axonly){ pty='n' } else { pty='p' }
    if(ptonly){ axs=FALSE } else { axs=TRUE }
    lh = TRUE 
    mark = t12marks[i]
    df = fullumap[[mark]]
    colnames(df)[colnames(df) == 'Project'] = 'project'
    df = df[order(df$GROUP, decreasing=TRUE),]
    xlim = range(df$X1)
    ylim = range(df$X2)
    # Imputed
    c = df[df$dataset == 'imp',]
    if (is.null(pchval)){ pch=19 } else { pch = pchvals[[pchval]][c[[pchval]]] }
    if (is.null(colval)){ col=c$COLOR } else { col = colvals[[colval]][c[[colval]]] }
    par(mar=c(1.5, 3.25, 0.25 + 1 * lh, 1))
    plot(c$X1, c$X2, col=col, pch=pch, yaxt='n', xaxt='n', axes=axs,
         ylab='', xlab='', cex=1, xlim=xlim, ylim=ylim, type=pty)
    if (!is.null(colval)){  
        print(colval)
        lcol = colvals[[colval]]
        colname = capitalize(colval)
        legend('bottomright', names(lcol), col=lcol, 
               pch=19, cex=1, inset=.01, box.col=NA, title=colname)
    }
    if (!ptonly){ 
        mtext(mark, side=2, cex=1.4, line=1.5)
        mtext('UMAP-1', side=1, cex=1, line=.25)
        mtext('UMAP-2', side=2, cex=1, line=0)
    }
    if (lh & (!ptonly)){ mtext('Imputed (all)', side=3, cex=1, col='black') }
    # Observed:
    c = df[df$dataset == 'obs',]
    obsid = as.character(c$id)
    if (is.null(pchval)){ pch=19 } else { pch = pchvals[[pchval]][c[[pchval]]] }
    if (is.null(colval)){ col=c$COLOR } else { col = colvals[[colval]][c[[colval]]] }
    par(mar=c(0.25, 3.25, 0.25 + 1 * lh, 0.5))
    plot(c$X1, c$X2, col=col, pch=pch, yaxt='n', xaxt='n', axes=axs,
         ylab='', xlab='', cex=.5, xlim=xlim, ylim=ylim, type=pty)
    if(lh & (!ptonly)) {mtext('Observed', side=3, cex=1, col='black')}
    # Only novel
    c = df[df$dataset == 'imp',]
    c = c[!(c$id %in% obsid),]
    if (is.null(pchval)){ pch=19 } else { pch = pchvals[[pchval]][c[[pchval]]] }
    if (is.null(colval)){ col=c$COLOR } else { col = colvals[[colval]][c[[colval]]] }
    par(mar=c(0.25, 2.25, 0.25 + 1 * lh, 1.5))
    plot(c$X1, c$X2, col=col, pch=pch, yaxt='n', xaxt='n', axes=axs,
         ylab='', xlab='', cex=.5, xlim=xlim, ylim=ylim, type=pty)
    if (lh & (!ptonly)){ mtext('Imputed (novel)', side=3, cex=1, col='black') }
    for (i in 2:NT12){ 
        mark = t12marks[i]
        df = fullumap[[mark]]
        colnames(df)[colnames(df) == 'Project'] = 'project'
        df = df[order(df$GROUP, decreasing=TRUE),]
        xlim = range(df$X1)
        ylim = range(df$X2)
        c = df[df$dataset == 'obs',]
        obsid = as.character(c$id)
        if (is.null(pchval)){ pch=19 } else { pch = pchvals[[pchval]][c[[pchval]]] }
        if (is.null(colval)){ col=c$COLOR } else { col = colvals[[colval]][c[[colval]]] }
        # Whether to put header or not:
        if ((i-1) %% 4 == 1){ lh = TRUE } else { lh = FALSE }
        # Observed:
        par(mar=c(0.25, 1.5, 0.25 + 1 * lh, 0.0))
        plot(c$X1, c$X2, col=col, pch=pch, yaxt='n', xaxt='n', axes=axs,
             ylab='', xlab='', cex=.5, xlim=xlim, ylim=ylim, type=pty)
        if (!ptonly) { mtext(mark, side=2, cex=1.4, line=0)}
        if(lh & (!ptonly)) {mtext('Observed', side=3, cex=1, col='black')}
        # Imputed
        c = df[df$dataset == 'imp',]
        if (is.null(pchval)){ pch=19 } else { pch = pchvals[[pchval]][c[[pchval]]] }
        if (is.null(colval)){ col=c$COLOR } else { col = colvals[[colval]][c[[colval]]] }
        par(mar=c(0.25, 0.0, 0.25 + 1 * lh, 0.0))
        plot(c$X1, c$X2, col=col, pch=pch, yaxt='n', xaxt='n', axes=axs,
             ylab='', xlab='', cex=.5, xlim=xlim, ylim=ylim, type=pty)
        if (lh & (!ptonly)){ mtext('Imputed (all)', side=3, cex=1, col='black') }
        # Only novel
        c = c[!(c$id %in% obsid),]
        if (is.null(pchval)){ pch=19 } else { pch = pchvals[[pchval]][c[[pchval]]] }
        if (is.null(colval)){ col=c$COLOR } else { col = colvals[[colval]][c[[colval]]] }
        par(mar=c(0.25, 0.0, 0.25 + 1 * lh, 0.25))
        plot(c$X1, c$X2, col=col, pch=pch, yaxt='n', xaxt='n', axes=axs,
             ylab='', xlab='', cex=.5, xlim=xlim, ylim=ylim, type=pty)
        if (lh & (!ptonly)){ mtext('Imputed (novel)', side=3, cex=1, col='black') }
    }
}


# png(paste0(imgpref,'umap_samerun_bothmarks_col_t12_largepanel_n',NNVAL, '.png'),res=300,units='in',width=17,height=6)
pdf(paste0(imgpref,'umap_samerun_bothmarks_col_t12_largepanel_n',NNVAL, '.pdf'),width=17,height=6)
plot.largepanel()
dev.off()

for (facet in names(pchvals)){
    png(paste0(imgpref,'umap_samerun_bothmarks_col_t12_largepanel_',facet,'_n',NNVAL, '.png'),res=300,units='in',width=17,height=6)
    plot.largepanel(facet)
    dev.off()
    pdf(paste0(imgpref,'umap_samerun_bothmarks_col_t12_largepanel_',facet,'_n',NNVAL, '.pdf'),width=17,height=6)
    plot.largepanel(facet)
    dev.off()
}

for (facet in names(colvals)){
    png(paste0(imgpref,'umap_samerun_bothmarks_col_t12_largepanel_col_',facet,'_n',NNVAL, '.png'),res=300,units='in',width=17,height=6)
    plot.largepanel(colval=facet)
    dev.off()
    pdf(paste0(imgpref,'umap_samerun_bothmarks_col_t12_largepanel_col_',facet,'_n',NNVAL, '.pdf'),width=17,height=6)
    plot.largepanel(colval=facet)
    dev.off()
}


# Plot points/axes separately:
pdf(paste0(imgpref,'umap_samerun_bothmarks_col_t12_largepanel_n',NNVAL, '_axonly.pdf'),width=17,height=6)
plot.largepanel(axonly=TRUE)
dev.off()

png(paste0(imgpref,'umap_samerun_bothmarks_col_t12_largepanel_n',NNVAL, '_ptonly.png'), res=500, units='in', width=17, height=6)
plot.largepanel(ptonly=TRUE)
dev.off()


# ---------------------------------------
# Add the read length and other metadata:
# ---------------------------------------
nam = 'all_submitted_released'
outfile <- paste0('Annotation/',nam,'_metadata.tsv')
metadata <- read.delim(outfile,sep="\t",header=T)
metadata$UID <- paste(metadata$Experiment.accession,
                      metadata$Biological.replicate.s.,
                      metadata$Technical.replicate,
                      sep=".")

fqmeta = metadata[metadata$File.format == 'fastq',]
# Year
fqmeta = unique(fqmeta[,c('Experiment.accession', 'Experiment.date.released', 'Lab',
                          'Run.type', 'Read.length', 'Platform', 'File.Status', 'Audit.NOT_COMPLIANT')])
names(fqmeta)[1] = 'Accession'
fqmeta = merge(accmap, fqmeta)
fqmeta$Year = sub("-.*","", fqmeta$Experiment.date.released)

rldf = aggregate(Read.length ~ Accession + id + Epitope, fqmeta, max)
col_fun = function(x, pal=colspec){
    palette = rev(pal)
    bin <- cut(x, seq(min(rldf$Read.length), max(rldf$Read.length), length.out=length(palette)), include.lowest=T) 
    palette[bin] }

# Plot the observed/imputed with metadata on the matched datasets:
png(paste0(imgpref,'umap_links_allmarks_impobs_n',NNVAL,'_readlength.png'), res=300, units='in', width=12, height=14)
layout(matrix(c(1:42),ncol=6,nrow=7, byrow=T) , heights=rep(2,7), widths=rep(2,6), TRUE)
sp=0
for (i in 1:(length(marks) -1)){ 
    mark = marks[i]
    udf = matchumap[[mark]]
    udf$COLOR = as.character(udf$COLOR)
    udf$COLOR[udf$COLOR == 'white'] <- 'lightgrey'
    udf$pch = 19
    udf$pch[udf$dataset == 'imp'] = 1
    sub.rldf = rldf[rldf$Epitope == mark,]
    udf = merge(udf, sub.rldf)
    # Separate out matched segments:
    count = aggregate(dataset ~ id, udf, function(x){length(unique(x))})
    keepct = as.character(count$id[count$dataset == 2])
    idf = merge(data.frame(id=keepct, dataset="imp"), udf[,c('id','dataset','X1','X2','COLOR')])
    colnames(idf)[3:4] = c('I1','I2')
    obdf = merge(data.frame(id=keepct, dataset="obs"), udf[,c('id','dataset','X1','X2','COLOR')])
    segdf = merge(obdf[,c('id','X1','X2','COLOR')], idf[,c('id','I1','I2')])
    c = udf
    c = c[order(c$GROUP, decreasing=TRUE),]
    # Print out ones with most len:
    segdf$len = sqrt((segdf$X1 - segdf$I1)^2 + (segdf$X2 - segdf$I2)^2)
    segdf = merge(segdf, meta)
    ssdf = segdf[segdf$len > 1,]
    xlim = c(min(c$X1), max(c$X1))
    ylim = c(min(c$X2), max(c$X2))
    # Make actual plot
    # c$COLOR = ifelse(c$dataset == 'imp','indianred','royalblue')
    # c$COLOR = col_fun(c$Read.length)
    par(mar=c(sp,sp,1.5,sp))
    plot(c$X1, c$X2, col=c$COLOR, pch=c$pch, ylab='UMAP-2', xlab='UMAP-1',
         cex=.75, xlim=xlim, ylim=ylim, axes=F)
    box(lwd=.5); mtext(paste0(mark))
    plot(c$X1, c$X2, col=ifelse(c$dataset == 'imp','indianred','royalblue'), pch=c$pch, ylab='UMAP-2', xlab='UMAP-1',
         cex=.75, xlim=xlim, ylim=ylim, axes=F)
    box(lwd=.5); mtext(paste0(mark))
    plot(c$X1, c$X2, col=col_fun(c$Read.length), pch=c$pch, ylab='UMAP-2', xlab='UMAP-1',
         cex=.75, xlim=xlim, ylim=ylim, axes=F)
    box(lwd=.5); mtext(paste0(mark))
}
plot(1,1, type='n', axes=F)
rls = sort(unique(rldf$Read.length))
legend('center', as.character(rls), pch=19, 
       col=col_fun(rls), pt.cex=2, cex=1.5, bty='n', ncol=3)
dev.off()


# ------------------------------------------------------------------------------------
# Show that modeling read length, lab, and other covariates 
# is barely accounting for variability when also accounting for lifestage, group, etc.
# ------------------------------------------------------------------------------------
# Prune fqmeta into a unique metadata table:
fqmeta = fqmeta[fqmeta$File.Status %in% c('released','submitted'),]
fqdf = aggregate(cbind(Read.length, Year) ~ Accession + Lab + Epitope + id + Run.type, fqmeta, max)
yrdf = aggregate(Year ~ Accession + Epitope + id, fqmeta, max)
lbdf = aggregate(Lab ~ Accession + Epitope + id, fqmeta, function(x){paste(sort(unique(x)), collapse=', ')})
rownames(rldf) = rldf$Accession
rownames(yrdf) = yrdf$Accession
rownames(lbdf) = lbdf$Accession


fdf = fqmeta[,c('Lab','Read.length','Year','Run.type')]

ggplot(fdf, aes(Read.length, Year, color=Lab)) + 
    geom_point() + theme_pubr()

ggplot(fdf, aes(Run.type, Read.length)) + 
    geom_boxplot() + theme_pubr()


fmat = cbind(rl=fqmeta$Read.length > 36, 
             yr=fqmeta$Year > 2013,
             pr=fqmeta$Project == 'Roadmap',
             rt=fqmeta$Run.type == 'single-ended')


# Attributes for datasets:
mdf = accmap[,c('id','Epitope','Accession')]
mdf$group = meta[mdf$id,'GROUP']
mdf$ls = meta[mdf$id,'lifestage']
mdf$type = meta[mdf$id,'type']
mdf$proj = meta[mdf$id,'Project']
mdf$rl = rldf[mdf$Accession,'Read.length']
mdf$yr = yrdf[mdf$Accession,'Year']
mdf$lb = lbdf[mdf$Accession,'Lab']


rmdf = aggregate(rl ~ group, mdf, mean)
rmdf = rmdf[order(rmdf$rl),]
mdf$group = factor(mdf$group, levels=rmdf$group)

ggplot(mdf, aes(group, rl, color=group)) + 
    geom_boxplot() + 
    geom_jitter() + 
    theme_pubr() + 
    scale_color_manual(values=colvals$group) + 
    theme(axis.text.x = element_text(angle=90, hjust=1), legend.position='none') 

mdf$rl.cut = cut(mdf$rl,breaks = c(0,36,50,75,90,101), right=FALSE)

# TODO: Move this all to a new script!
fmat = cbind(rl=mdf$rl> 36, 
             yr=mdf$yr > 2013,
             pr=mdf$rl == 'Roadmap',
             # rt=mdf$rt == 'single-ended',
             tp=mdf$type == 'tissue')
dist(t(fmat))

# TODO: Enrichments by lab - all 




# Gather all distances as a table
aodf = c()
aidf = c()
for (i in 1:(length(marks) -1)){
    # Load both observed and imputed datasets:
    mark = marks[i]
    print(mark)
    omat = obsll[[mark]]
    imat = ll[[mark]]
    # Restrict to chosen cells:
    cn = rownames(omat)[rownames(omat) %in% cellorder]
    omat = omat[cn, cn]
    icn = rownames(imat)[rownames(imat) %in% cellorder]
    imat = imat[icn, icn]
    # Also calculate the same exact restricted matrix:
    bcn = icn[icn %in% cn]
    names(bcn) = NULL
    bomat = omat[bcn, bcn]
    bimat = imat[bcn, bcn]
    bodf = data.frame(bomat * upper.tri(bomat))
    bodf$id = rownames(bodf)
    bodf = gather(bodf, id2, dist, -id)
    bodf$Epitope = mark
    bidf = data.frame(bimat * upper.tri(bimat))
    bidf$id = rownames(bidf)
    bidf = gather(bidf, id2, dist, -id)
    bidf$Epitope = mark
    # Add: 
    aodf = rbind(aodf, bodf)
    aidf = rbind(aidf, bidf)
}
aodf = aodf[aodf$id != aodf$id2,]
aidf = aidf[aidf$id != aidf$id2,]
aodf = aodf[aidf$dist != 0,] # Remove double counting.
aidf = aidf[aidf$dist != 0,]
aodf$type = 'Observed'
aidf$type = 'Imputed'

# Need to model distance by is same of each of x attributes.
accmap1 = accmap[,c('id','Epitope','Accession')]
accmap2 = accmap1
names(accmap2) = c('id2', 'Epitope','Accession2')
aodf = merge(merge(aodf, accmap1), accmap2)
aidf = merge(merge(aidf, accmap1), accmap2)

# Add observations:
aodf$same.group = meta[aodf$id,'GROUP'] == meta[aodf$id2, 'GROUP']
aodf$same.ls = meta[aodf$id,'lifestage'] == meta[aodf$id2, 'lifestage']
aodf$same.type = meta[aodf$id,'type'] == meta[aodf$id2, 'type']
aodf$same.proj = meta[aodf$id,'Project'] == meta[aodf$id2, 'Project']
aodf$same.rl = rldf[aodf$Accession,'Read.length'] == rldf[aodf$Accession2, 'Read.length']
aodf$same.yr = yrdf[aodf$Accession,'Year'] == yrdf[aodf$Accession2, 'Year']
aodf$same.lb = lbdf[aodf$Accession,'Lab'] == lbdf[aodf$Accession2, 'Lab']

# Add observations:
aidf$same.group = meta[aidf$id,'GROUP'] == meta[aidf$id2, 'GROUP']
aidf$same.ls = meta[aidf$id,'lifestage'] == meta[aidf$id2, 'lifestage']
aidf$same.type = meta[aidf$id,'type'] == meta[aidf$id2, 'type']
aidf$same.proj = meta[aidf$id,'Project'] == meta[aidf$id2, 'Project']
aidf$same.rl = rldf[aidf$Accession,'Read.length'] == rldf[aidf$Accession2, 'Read.length']
aidf$same.yr = yrdf[aidf$Accession,'Year'] == yrdf[aidf$Accession2, 'Year']
aidf$same.lb = lbdf[aidf$Accession,'Lab'] == lbdf[aidf$Accession2, 'Lab']

# i0 = lm(dist ~ same.group + same.ls + same.type + same.proj + same.yr + Epitope, aidf)
# i1 = lm(dist ~ same.group + same.ls + same.type + same.proj + same.rl + same.yr + Epitope, aidf)
# summary(i1)
# anova(i0,i1)

# o0 = lm(dist ~ same.group + same.ls + same.type + same.proj + same.yr + Epitope, aodf)
# o1 = lm(dist ~ same.group + same.ls + same.type + same.proj + same.rl + same.yr + Epitope, aodf)
# summary(o1)
# anova(o0, o1)

mark = 'H3K27ac'
o0 = lm(dist ~ same.group + same.ls + same.type + same.proj + same.yr, aodf[aodf$Epitope == mark, ])
o1 = lm(dist ~ same.group + same.ls + same.type + same.proj + same.rl + same.yr, aodf[aodf$Epitope == mark, ])
summary(o1)
anova(o0, o1)

i0 = lm(dist ~ same.group + same.ls + same.type + same.proj + same.yr, aidf[aidf$Epitope == mark, ])
i1 = lm(dist ~ same.group + same.ls + same.type + same.proj + same.rl + same.yr, aidf[aidf$Epitope == mark, ])
summary(i1)
anova(i0, i1)


mark = 'H3K27ac'
coeffdf = c()
for (i in 1:(length(marks) - 1)){
    mark = marks[i]
    print(mark) 
    bdf = rbind(aodf[aodf$Epitope == mark, ], aidf[aidf$Epitope == mark,])
    # i0 = lm(dist ~ same.group + same.ls + same.type + same.proj + same.yr, bdf)
    i1 = lm(dist ~ type*(same.group + same.ls + same.type + same.proj + same.rl + same.yr + same.lb), bdf)
    cfdf = data.frame(coefficients(summary(i1)))
    colnames(cfdf) = c('Est', 'SE','tvalue', 'pvalue')
    cfdf$term = rownames(cfdf)
    cfdf$mark = mark
    coeffdf = rbind(coeffdf, cfdf)
    # av = anova(i0, i1)
}

coeffdf$log10p = -log10(coeffdf$pvalue)
ewide = spread(coeffdf[, c('Est','mark','term')], mark, Est)
pwide = spread(coeffdf[, c('log10p','mark','term')], mark, log10p)
emat = as.matrix(ewide[,-1])
rownames(emat) = ewide$term
pmat = as.matrix(pwide[,-1])
rownames(pmat) = pwide$term

emat = emat[grep('same', rownames(emat)),]
pmat = pmat[grep('same', rownames(pmat)),]
cov.map = c('Year','Sample Type','Read Len','Project','Lifestage','Lab','Bio. group')
names(cov.map) = c('yr','type','rl','proj','ls','lb','group')

png(paste0(imgpref,'covars_regression_allmarks_readlength.png'), res=300, units='in', width=5, height=5)
mx = max(abs(emat), na.rm=T)
par(mar=c(6,6,0,0))
image(t(emat), col=rev(colrb), zlim=c(-mx, mx), axes=F)
box()
ind = which(pmat > -log10(.05),arr.ind=TRUE)
vals = pmat[ind]
sig = ifelse(vals > -log10(.01), ifelse(vals > -log10(0.001), '***', '**'), '*')
xat = seq(0,1, length.out=ncol(emat))
yat = seq(0,1, length.out=nrow(emat))
covstr = sub("TRUE", "", rownames(emat))
covstr = sub("^type", "", covstr)
oind = grep("Observed", rownames(emat))
# covstr = sub("Observed:", "Obs: ", covstr)
covstr = sub("Observed:", "", covstr)
# covstr = sub("^same", "Imp: same", covstr)
covstr = sub("same.", "", covstr)
covstr = cov.map[covstr]
abline(h=0.5, xpd=TRUE)
text(x=xat[ind[,2]], y=yat[ind[,1]], 
     sig, xpd=TRUE, adj =.5)
text(x=parpos(1, .02), y=yat, 
     covstr, xpd=TRUE, adj =1, 
     col=ifelse(1:nrow(emat) %in% oind, 'royalblue','indianred'))
text(y=parpos(2, .02), x=xat, 
     colnames(emat), xpd=TRUE, adj =1, srt=90)
dev.off()


# Alternatively, z-score the per-mark distances (in observed and imputed)
bdf = rbind(aodf, aidf)
mbdf = aggregate(dist ~ Epitope + type, bdf, mean)
sbdf = aggregate(dist ~ Epitope + type, bdf, sd)
names(mbdf)[3] = 'mean'
names(sbdf)[3] = 'sd'
bdf = merge(bdf, merge(mbdf, sbdf))

bdf$zscore = (bdf$dist - bdf$mean) / bdf$sd
i1 = lm(dist ~ type*(same.group + same.ls + same.type + same.proj + same.rl + same.yr + same.lb), bdf)
cfdf = data.frame(coefficients(summary(i1)))
colnames(cfdf) = c('Est', 'SE','tvalue', 'pvalue')
cfdf$term = rownames(cfdf)
cfdf = cfdf[grep('same', rownames(cfdf)),]

png(paste0(imgpref,'covars_regression_allmarks_readlength_zscored.png'), res=300, units='in', width=2, height=4)
mx = max(abs(cfdf$Est), na.rm=T)
sp =0.1
par(mar=c(sp,8,sp,sp))
image(t(cfdf$Est), col=rev(colrb), zlim=c(-mx, mx), axes=F)
box()
vals = -log10(cfdf$pvalue)
xat = parpos(1,-.5)
yat = seq(0,1, length.out=nrow(cfdf))
oind = grep("Observed", rownames(cfdf))
iind = 1:nrow(cfdf)
iind = iind[!(iind %in% oind)]
covstr = sub("TRUE", "", rownames(cfdf))
covstr = sub("^type", "", covstr)
covstr = sub("Observed:", "", covstr)
covstr = sub("same.", "", covstr)
covstr = cov.map[covstr]
abline(h=0.5)
text(x=xat, y=yat, round(vals,0), xpd=TRUE, adj =.5)
text(x=parpos(1, .1), y=yat, 
     covstr, xpd=TRUE, adj =1, 
     col=ifelse(1:nrow(emat) %in% oind, 'royalblue','indianred'))
text(x=parpos(1, 3.5), y=mean(yat[oind]), "Observed", xpd=TRUE, adj =.5, srt=90)
text(x=parpos(1, 3.5), y=mean(yat[iind]), "Imputed", xpd=TRUE, adj =.5, srt=90)
ypad=0.01
segments(x0=parpos(1,3), y0=yat[c(min(iind), min(oind))] + ypad,
         x1=parpos(1,3), y1=yat[c(max(iind), max(oind))] - ypad, xpd=TRUE)
dev.off()


# ----------------------------------
# Plot reduced panels for figure 2b:
# ----------------------------------
redmarks = c("H3K27ac", "H3K27me3",
             "H3K4me3","H3K4me1","H3K36me3", 
             "DNase-seq", "H3K9me3", "H3K9ac")

plot.redpanel = function(pchval=NULL, axonly=FALSE, ptonly=FALSE){
    lmat = matrix(c(1,1, 2,2, 3,4,5, 1,1, 2,2, 6,7,8), nrow=2, byrow=T)
    layout(lmat, heights=c(1.1,1.1), widths=rep(1, 7), TRUE)
    # i = 1
    if(axonly){ pty='n' } else { pty='p' }
    if(ptonly){ axs=FALSE } else { axs=TRUE }
    if (is.null(pchval)){ pch=19 } else { pch = pchvals[[pchval]][c[[pchval]]] }
    for (i in 1:length(redmarks)){ 
        mark = redmarks[i]
        df = fullumap[[mark]]
        df = df[order(df$GROUP, decreasing=TRUE),]
        xlim = range(df$X1)
        ylim = range(df$X2)
        c = df[df$dataset == 'imp',]
        if (is.null(pchval)){ pch=19 } else { pch = pchvals[[pchval]][c[[pchval]]] }
        par(mar=c(0, 0.0, 0 + 1 * lh, 0.0))
        ptcex = ifelse(i <= 2, 0.65, 0.4)
        plot(c$X1, c$X2, col=c$COLOR, pch=pch, yaxt='n', xaxt='n', axes=axs,
             ylab='', xlab='', cex=ptcex, xlim=xlim, ylim=ylim, type=pty)
        if (lh & (!ptonly)){ mtext(mark, side=3, cex=1, col='black') }
        # Calculate/plot arrows:
        cm = apply(c[, c('X1', 'X2')], 2, mean)
        cmdf = c()
        # covariates = colnames(metamat)[1:3]
        covariates = colnames(metamat)[c(1,3)]
        tolab = rbind(data.frame(cov='type', facet=c('tissue', 'cell line', 'primary cell')), 
                      data.frame(cov='lifestage', facet=c('embryonic', 'adult')), 
                      data.frame(cov='group', facet=c('Blood & T-cell', 'Brain','Muscle')))
        colnames(c)[6] = 'group'
        # for (cov in covariates){
        #     for (facet in unique(c[[cov]])){
        for (j in 1:nrow(tolab)){
            cov = as.character(tolab[j, 'cov'])
            facet = as.character(tolab[j, 'facet'])
            cm2 = apply(c[c[[cov]]== facet, c('X1', 'X2')], 2, mean)
            cmdf = rbind(cmdf, data.frame(x0=cm[1], y0=cm[2], 
                                          x1=cm2[1], y1=cm2[2], 
                                          cov=cov, facet=facet, col=colvals[[cov]][facet]))
        }
        # }
        # }
        cmdf$col = as.character(cmdf$col)
        cmdf$facet = sapply(as.character(cmdf$facet), capitalize)
        cmdf = cmdf[cmdf$col != 'white',]
        require(sfsmisc)
        if (axonly){
            p.arrows2(x1=cmdf$x0, y1=cmdf$y0, x2=cmdf$x1, y2=cmdf$y1, 
                      col=cmdf$col, fill=cmdf$col, size=ifelse(i <=2, .7, .5), width=1.5, lwd=1, border=NULL)
            lbcex=ifelse(i <=2, 0.75, 0.4)
            rdf = general_repel_text(x=cmdf$x1, y=cmdf$y1, labels=cmdf$facet,
                                     xlim=par()$usr[1:2], ylim=par()$usr[3:4],
                                     hjust=.5, vjust=.5, seed=1, max.iter=5000,
                                     cex=lbcex, pt.cex=.25)
            text(x=rdf$x, y=rdf$y, labels=rdf$lab,
                 srt=0, adj=0, xpd=TRUE, cex=lbcex, col=cmdf$col)
        }
    }
}


# Plot points/axes separately:
pdf(paste0(imgpref,'umap_samerun_bothmarks_col_t12_redpanel_n',NNVAL, '_axonly.pdf'),width=8,height=3)
plot.redpanel(axonly=TRUE)
dev.off()

png(paste0(imgpref,'umap_samerun_bothmarks_col_t12_redpanel_n',NNVAL, '_ptonly.png'), res=500, units='in', width=8, height=3)
plot.redpanel(ptonly=TRUE)
dev.off()






# --------------------------------------
# Plot UMAP with imputed-observed links:
# --------------------------------------
for (i in 1:(length(marks) -1)){
    mark = marks[i]
    mat = full.ll[[mark]]
    # Restrict to chosen cells:
    samp = t(sapply(rownames(mat),
                    function(x){strsplit(x, "_")[[1]][1]}))
    cn = rownames(mat)[samp %in% cellorder]
    mat = mat[cn, cn]
    # Use previous UMAP:
    # nn = 100
    # mdist = 0.3
    # u = umap(as.dist(mat), n_neighbors=nn, min_dist=mdist, verbose=F, repulsion_strength=.25)
    udf = fullumap[[mark]]

    # Annotation:
    nam = colnames(mat)
    mdf = t(sapply(nam, function(x){strsplit(x, "_")[[1]]}))
    udf$COLOR = as.character(udf$COLOR)
    udf$COLOR[udf$COLOR == 'white'] <- 'lightgrey'
    udf$pch = 19
    udf$pch[udf$dataset == 'imp'] = 1

    # Separate out matched segments:
    count = aggregate(dataset ~ id, udf, function(x){length(unique(x))})
    keepct = as.character(count$id[count$dataset == 2])
    idf = merge(data.frame(id=keepct, dataset="imp"), udf[,c('id','dataset','X1','X2','COLOR')])
    colnames(idf)[3:4] = c('I1','I2')
    odf = merge(data.frame(id=keepct, dataset="obs"), udf[,c('id','dataset','X1','X2','COLOR')])
    segdf = merge(odf[,c('id','X1','X2','COLOR')], idf[,c('id','I1','I2')])
    c = udf
    c = c[order(c$GROUP, decreasing=TRUE),]

    # TODO: PRINT OUT ONES WITH MOST LEN:
    segdf$len = sqrt((segdf$X1 - segdf$I1)^2 + (segdf$X2 - segdf$I2)^2)
    segdf = merge(segdf, meta)
    ssdf = segdf[segdf$len > 1,]

    xlim = c(min(c$X1), max(c$X1))
    ylim = c(min(c$X2), max(c$X2))

    png(paste0(imgpref,'umap_links_', mark, '_n',NNVAL,'.png'),res=300,units='in',width=10,height=8)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    # Make actual plot
    plot(c$X1, c$X2, col=c$COLOR, pch=c$pch, ylab='UMAP-2', xlab='UMAP-1',cex=.5, xlim=xlim, ylim=ylim)
    segments(segdf$X1,segdf$X2, segdf$I1, segdf$I2, col=segdf$COLOR)
    mtext(paste0('UMAP on ', metric,' of ', mark, '\n Imputed (open) and Observed (closed)'), side=3, line=1)
    text(ssdf$X1, ssdf$X2, labels=paste0(ssdf$id, "(", ssdf$GROUP,")\n", ssdf$infoline), col=ssdf$COLOR, cex=.5)
    # Add the legend:
    upViewport()
    draw(col.list[['group']], x = circle_size, just = "left")
    dev.off()

    # Highlight only the linked points:
    kid = which(c$id %in% keepct)
    png(paste0(imgpref,'umap_links_ol_', mark, '_n', NNVAL, '.png'),res=300,units='in',width=10,height=8)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    # Make actual plot
    plot(c$X1[-kid], c$X2[-kid], col='lightgray', pch=c$pch[-kid], ylab='UMAP-2', xlab='UMAP-1', 
         cex=.5, xlim=xlim, ylim=ylim)
    points(c$X1[kid], c$X2[kid], col=c$COLOR[kid], pch=c$pch[kid], ylab='UMAP-2', xlab='UMAP-1',cex=.5)
    segments(segdf$X1,segdf$X2, segdf$I1, segdf$I2, col=segdf$COLOR)
    mtext(paste0('UMAP on ', metric, ' of ', mark, '\n Imputed (open) and Observed (closed)'), side=3, line=1)
    text(ssdf$X1, ssdf$X2, labels=paste0(ssdf$id, "(", ssdf$GROUP,")\n", ssdf$infoline), col=ssdf$COLOR, cex=.5)
    # Add the legend:
    upViewport()
    draw(col.list[['group']], x = circle_size, just = "left")
    dev.off()
}

