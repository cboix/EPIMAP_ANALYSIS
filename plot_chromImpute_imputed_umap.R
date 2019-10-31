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


# UMAP from full data:
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

