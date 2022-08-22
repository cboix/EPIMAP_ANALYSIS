#!/usr/bin/R
# -----------------------------------------------------
# Load modules from CHMM/other clustering
# + get three different orderings
# -----------------------------------------------------
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
library(proxy)  # Needed for ejaccard dist
library(cba)
library(viridis)
library(tidyr)
library(dplyr)
library(scales)

# Arguments:
# filepref = 'cls_300'
# tagline = 'ChromHMM Enhancers'
# imgdir = paste0(img, "clusters/") 
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need clusters filename")
} else {        
    filepref = args[1]
    tagline = args[2]
    if (length(args) > 2){ 
        imgdir = args[3] 
    }
}

# Prefix:
cmd = paste('mkdir -p', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'clusters_')

# Relevant files:
centersfile = paste0(filepref, '_centers.tsv.gz')
countsfile = paste0(filepref, '_counts.tsv')
bpfile = paste0(filepref, '_nbp.tsv')
namesfile = paste0(filepref, '_names.tsv')

# --------------
# Load datasets:
# --------------
centers = as.matrix(read.delim(gzfile(centersfile),header=F))
colnames(centers) = paste0('c', 0:(ncol(centers) - 1))
rownames(centers) = scan(namesfile,'c')
print(paste0("Loaded file with shape: ",
             paste0(dim(centers), collapse=", ")))

# Reduce to epigenomes. 
# Split and merge names if necessary. 
# TODO: Also create plots with state?
# TODO: Plot each chmm element with appropriate colors?
filtcutoff = 0.01
centdf = data.frame(centers, nam=rownames(centers))
clong = gather(centdf, cls, value, -nam)
clong = filter(clong, clong$value > filtcutoff)
rn = rownames(centers)
rnsplit = sapply(rn, function(x) strsplit(x,"_")[[1]][1])
rnsplit2 = sapply(rn, function(x) strsplit(x,"_")[[1]][2])
reducenames = data.frame(nam=rn, id=rnsplit, state=rnsplit2)
clong = merge(clong, reducenames)

# Reduce:
cred = aggregate(value ~ id + cls, clong, max)
cwide = spread(cred, cls, value, fill=0)
cmat = as.matrix(cwide[,-1])
rownames(cmat) = cwide$id

# Metadata:
ctdf = read.delim(countsfile, header=F, stringsAsFactors=F)
bpdf = read.delim(bpfile, header=F, stringsAsFactors=F)
names(bpdf) = c('cls','nbp')
names(ctdf) = c('cls','nregion')
clsdf = merge(ctdf, bpdf)
clsdf$avg.bp = clsdf$nbp / clsdf$nregion  # Average of 200

# For plotting counts:
ct = clsdf$nregion
names(ct) = paste0('c', clsdf$cls)

# As lists:
counts = list()
centlist = list()
counts[[tagline]] = ct
centlist[[tagline]] = cmat

# ---------------------------------------------------
# Hypergeometric testing for enrichment for clusters:
# ---------------------------------------------------
# Run hypergeometric for each comb, make table:
# TODO: Add Enrichment against metadata 
# - Aka which clusters are more adult vs. which top (1-2) CT
# For each cluster center, evaluate against each metadata column
# TODO: Add also the availability of different marks as enrichment factor
# ---------------------------------------------------
# Generate for all:
threshold = 0.1
covariates = c('lifestage','sex','type','Project','GROUP')
# For run.hyper give vector with hits, draws, white balls, total balls:
run.hyper <- function(y){phyper(q=y[1]-1, m=y[3], n=y[4] - y[3], k=y[2], lower.tail=FALSE)} 
run.fisher <- function(y){
    table <-  matrix(as.numeric(c(y[1], y[2], y[3], y[4])), ncol = 2, byrow = TRUE)
    if(any(is.na(table))) p <- "error" else p <- fisher.test(table, alternative="two.sided")$p.value
    return(p)
}

hdf = filter(cred, value > threshold)
hdf = merge(hdf, meta[,c('id',covariates)])
# For each cluster (relabel) - test
# keptid = levels(reducenames$id)
keptid = unique(reducenames$id)
submeta = meta[keptid,]

hgdf = NULL
subhgdf = NULL
redhgdf = NULL
enrbreaks = c(0)
for (cov in covariates){
    lvls = unique(hdf[[cov]])
    # lvls = levels(hdf[[cov]])
    cat(paste0('COVARIATE: ', cov, 
               '\nLEVELS: ', paste0(lvls, collapse=", "), "\n"))
    enrbreaks = c(enrbreaks, enrbreaks[length(enrbreaks)] + length(lvls))
    # Do multiple fishers or a hypergeom?
    # TEST: For single cluster, which covariate levels are more enriched - 1 vs. 0 
    # For each level, test:
    for (lv in lvls){
        draws = table(hdf[, 'cls'])
        subhits = table(hdf[hdf[[cov]] == lv, 'cls'])
        hits = draws * 0
        hits[names(subhits)] = subhits
        df = cbind(q=hits, draw=draws,
                   m=nrow(submeta[submeta[[cov]] == lv,]),
                   N=nrow(submeta))
        pout <- apply(df, 1, run.hyper)
        df = data.frame(df, p.value=pout, cls=rownames(df), covariate=cov, level=lv)
        rownames(df) = NULL
        if (is.null(hgdf)){hgdf = df} else {hgdf = rbind(hgdf, df)}
    }
}

# Turn into wide matrix and plot:
hgdf$log10p = -log10(hgdf$p.value)
hgdf$log10p[hgdf$p.value < 10^-12] <- 12  # Cap inf
hgdf$cl = paste0(hgdf$covariate, '\n', hgdf$level)
hgwide = spread(hgdf[, c('cls','covariate','cl','log10p')], cls, log10p)
hgmat = as.matrix(hgwide[,-c(1,2)])
CUTOFF=5
hgmat[hgmat > CUTOFF] <- CUTOFF
hgmat[hgmat < 2] <- 0
rownames(hgmat) = hgwide$cl

# Reduce the plot to top 2-3 cells + top of others.
# TODO Add histones
enrichmat = matrix("unknown/mixed", nrow=ncol(cmat), ncol=ncol(metamat),
                   dimnames = list(colnames(cmat), colnames(metamat)))
for (cov in covariates){
    # Min p.value per cluster:
    sdf = filter(hgdf, covariate == cov)
    mdf = aggregate(p.value ~ cls, sdf, min)
    mdf = merge(mdf, sdf)
    # Remove all log10p lower than 2
    mdf = filter(mdf, log10p > 2)
    # Reduce to one per (tie breaking by lvl)
    mdf = aggregate(level ~ cls,  mdf, function(x){head(x,1)})
    covname = tolower(cov)
    enrichmat[as.character(mdf$cls), covname] = as.character(mdf$level)
}

# ---------------------------
# Plot clusters side by side:
# ---------------------------
# TODO: Make these functions standalone
plot.centers = function(centers, set, cellorder, counts, calc.only=FALSE, title=TRUE,
                        subset=FALSE, cls=acut.nam, cls.ord=NULL, palette=col1, ablwd=1){
    mat = centers[[set]]
    ct = counts[[set]]
    if (subset) {
        cidx = keep.cls[[set]]
        mat = mat[,cidx]
        ct = ct[cidx]
    }
    idx = as.numeric(which(ct != 0))
    mat = mat[,idx]
    if (is.null(cls.ord)){
        subid = cellorder[cellorder %in% rownames(mat)]
        tmp = matrix(NA, nrow=length(cellorder), ncol=ncol(mat))
        rownames(tmp) = cellorder
        # Reorder:
        ll = diag.mat(mat[subid,])
        tmp[subid,]= ll[[1]]
        # Use cto to get breaks:
        cto = ll[[3]]
        vcut =  c(cto[cto ==0], cls[subid][cto])
        vbreaks = calc.breaks.acut(vcut)
    } else {
        subid = cellorder[cellorder %in% rownames(mat)]
        subclsord = cls.ord[cls.ord %in% colnames(mat)]
        tmp = matrix(NA, nrow=length(cellorder), ncol=length(subclsord),
                     dimnames=list(cellorder, subclsord))
        # Reorder:
        tmp[subid,subclsord] = mat[subid,subclsord]
        ll = list(subclsord, subclsord)
        vcut = round(1:ncol(tmp) / 20) 
        vbreaks = calc.breaks.acut(vcut)
    }
    if (!calc.only){
        image(t(tmp), axes=F, col=palette, zlim=c(0,1), useRaster=TRUE)
        abline(h=par()$usr[3:4],lty=1,lw=0.5)
        abline(v=par()$usr[1:2],lty=1,lw=0.5)
        abline(v=vbreaks,lty='dotted',lw=ablwd, col='darkgrey')
        if (title){ mtext(set, side=3, cex=1.3) }
    }
    return(list(ll[[2]], vbreaks, vcut))
}

plot.counts = function(counts, set, ordll, ylim=NULL, subset=FALSE, cex.axis=0.75){
    ct = counts[[set]]
    ord = ordll[[1]]
    vbreaks = ordll[[2]] # * length(ct)
    # Using 
    if (subset) { ct = ct[keep.cls[[set]]] }
    ct = ct[ct != 0]
    subord = ord[ord %in% names(ct)]
    ct = ct[subord]
    barplot(-ct, ylim=ylim, col='darkgrey', border=NA, yaxt='n',xaxt='n')
    at = seq(0,-80, length.out=5) * 1000
    lbs = comma_format()(seq(0, 80, length.out=5) * 1000)
    axis(2, lwd=.5, las=1, at=at, labels=lbs, cex.axis=cex.axis)
    # abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
}


plot.countprop = function(counts, set, ordll, ylim=NULL, 
                          subset=FALSE, cex.axis=0.75, cols=NULL){
    ct = counts[[set]]
    ord = ordll[[1]]
    vbreaks = ordll[[2]] # * length(ct)
    # Reorder:
    if (subset) { ct = ct[keep.cls[[set]]] }
    ct = ct[ct != 0]
    subord = ord[ord %in% names(ct)]
    ct = ct[subord]
    # Create cuts:
    tot = sum(ct)
    x1 = c(0,cumsum(ct)) / tot
    x2 = seq(0, 1, length.out = length(x1))
    # Plot:
    par(yaxs='i')
    par(xaxs='i')
    plot(0, type='n', xlim=c(0,1), ylim=c(0,3), axes=FALSE)
    # Polygons:
    linecols = rep('black', length(x1))
    if (!is.null(cols)){
        for (i in 1:(length(x1)-1)){
            polygon(x= c(x1[i], x1[i+1], x1[i+1], x2[i+1], x2[i+1], x2[i], x2[i], x1[i]),
                    y=c(0, 0, 1, 2.25, 3, 3, 2.25, 1), col=cols[i], border=NA)
            if (cols[i] == "#000000") { 
                linecols[c(i, i+1)] = 'grey25'
            }
        }
    }
    # Segments:
    segments(x0=x1, y0=0, x1=x1, y1=1, lwd=.25, col=linecols)
    segments(x0=x1, y0=1, x1=x2, y1=2.25, lwd=.25, col=linecols)
    segments(x0=x2, y0=2.25, x1=x2, y1=3, lwd=.25, col=linecols)
    box(lwd=.25)
    sq = seq(0, tot, 1e5)
    at = seq(0,sq[length(sq)] / tot, length.out=length(sq))
    lbs = comma_format()(sq)
    axis(1, lwd=.25, las=1, at=at, labels=rep('', length(at)))
    text(x=at, y=par()$usr[3] - 0.5 * diff(par()$usr[3:4]),
         labels=lbs, cex=cex.axis, xpd=TRUE, srt=0)
}


plot_clusters <- function(rnorder, acutnamed, lablist, plot.rect=FALSE) {
    layout(matrix(c(1:9),3,3), heights=c(8,.3,1), widths=c(.9,8,.5), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    dist.order = list()
    # Metadata matrix:
    par(mar=c(0.25, 6, 2, 0))
    meta.image(metamat[rnorder,5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    text(y=lablist[[1]],
         x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
         labels=lablist[[2]], srt=0, adj=1, xpd=TRUE,cex=.8, col=lablist[[3]])
    abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
    # Add labels
    par(mar=c(0.25, 6, 0, 0))
    image(t(as.matrix(metamat) == ""), col='white', axes=F, useRaster=TRUE)
    metaclass = sapply(rev(colnames(metamat)), capitalize)
    text(x=seq(0,1, length.out=ncol(metamat)),
         y=par()$usr[4]-0.01*(par()$usr[4]-par()$usr[3]), 
         labels=metaclass, srt=90, adj=1, xpd=TRUE, cex=.7)
    plot.new()
    # Plot clusters and counts
    set = tagline
    par(mar=c(0.25, 0.25, 2, 0.25))
    dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed)
    abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
    if (plot.rect){
        # Add rectangles to centers:
        rll = calc.breaks.rect(hcls=acutnamed, vcls=dist.order[[set]][[3]], colset)
        rectdf = rll[[1]]
        vccols = rll[[2]]
        # Add grey - ubq rectangle:
        rectdf = rbind(c(x1=par()$usr[1], x2=rectdf$x1[1],
                         y1=par()$usr[3], y2=par()$usr[4]), rectdf)
        vccols = c('grey',vccols)
        rect(xleft=rectdf$x1, xright=rectdf$x2,
             ybottom=rectdf$y1, ytop=rectdf$y2,
             border=vccols, lwd=1)
        par(xpd=NA)
        rect(xleft=rectdf$x1, xright=rectdf$x2, ybottom=par()$usr[3], 
             ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
             border='white', col=vccols, lwd=.25)
        par(xpd=FALSE)
    }
    par(mar=c(0.25, 0.25, 0, 0.25))
    clsord = dist.order[[set]][[1]]
    meta.image(enrichmat[clsord,1:5], colvals=colvals, cex=0, horiz=F, useRaster=TRUE)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    vbreaks = dist.order[[set]][[2]]
    abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.counts(counts, set, dist.order[[set]])
    # Availability:
    par(mar=c(.25, 0, 2, 0.25))
    avail = as.matrix(wm[main.marks, rnorder])
    image(avail, axes=F, col=c('white', 'darkgrey'), useRaster=TRUE)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    par(mar=c(.25, 0, 0, 0.25))
    image(avail, axes=F, col='white', useRaster=TRUE)
    text(x=seq(0,1, length.out=length(main.marks)),
         y=par()$usr[4]-0.01*(par()$usr[4]-par()$usr[3]), 
         labels=main.marks, srt=90, adj=1, xpd=TRUE,cex=.65)
    return(dist.order[[set]])
}

plot_clusters_vsmall <- function(rnorder, acutnamed, lablist, plot.rect=FALSE, small=FALSE, tiny=FALSE, count.barplot=TRUE) {
    if (small){
        layout(matrix(c(1:6),3,2), heights=c(5,.4,.6), widths=c(1.5,8), TRUE)
    } else if (tiny) {
        # layout(matrix(c(1:6),3,2), heights=c(5,.4,.7), widths=c(1.5,8), TRUE)
        layout(matrix(c(1:6),3,2), heights=c(5,.4,.6), widths=c(2,8), FALSE)
        small=TRUE
    } else {
        layout(matrix(c(1:6),3,2), heights=c(8,.5,.8), widths=c(1.7,8), TRUE)
    }
    par(yaxs="i")
    par(xaxs="i")
    dist.order = list()
    # Metadata matrix:
    par(mar=c(0.25, 7 + 2 * tiny, .25, 0))
    meta.image(metamat[rnorder,5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    # TODO: make small text repel function:
    box.pad = 0.018 + 0.002 * small  + 0.003 * tiny
    xx = space.1d(lablist[[1]], box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
    xx = space.1d(xx, box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
    x0 = par()$usr[1]-0.1*(diff(par()$usr[1:2]))
    x1 = par()$usr[1]-0.15*(diff(par()$usr[1:2]))
    x2 = par()$usr[1]-0.3*(diff(par()$usr[1:2]))
    x3 = par()$usr[1]-0.35*(diff(par()$usr[1:2]))
    text(y=xx, x=x3, labels=lablist[[2]],
         srt=0, adj=1, xpd=TRUE, 
         cex=1 - 0.2 * small - 0.1 * tiny, col=lablist[[3]])
    par(xpd=TRUE)
    segments(x0=x2, y0=xx, x1=x1, y1=lablist[[1]], col=lablist[[3]])
    segments(x0=x0, y0=lablist[[1]], x1=x1, y1=lablist[[1]], col=lablist[[3]])
    par(xpd=FALSE)
    abline(h=dist.breaks,lty='dotted', lw=.5, col='darkgrey')
    # Add labels
    par(mar=c(0.25, 7 + 2 * tiny, 0, 0))
    image(t(as.matrix(metamat) == ""), col='white', axes=F, useRaster=TRUE)
    metaclass = sapply(rev(colnames(metamat)), capitalize)
    text(x=seq(0,1, length.out=ncol(metamat)),
         y=par()$usr[4]-0.01*(par()$usr[4]-par()$usr[3]), 
         labels=metaclass, srt=90, adj=1, xpd=TRUE, cex=.65 - .1 * tiny)
    # Blank
    plot.new()
    # Plot clusters and counts
    set = tagline
    par(mar=c(0.25, 0.25, 0.25, 0.25))
    dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, title=FALSE, ablwd=.5)
    abline(h=dist.breaks,lty='dotted', lw=.5, col='darkgrey')
    if (plot.rect){
        # Add rectangles to centers:
        rll = calc.breaks.rect(hcls=acutnamed, vcls=dist.order[[set]][[3]], colset)
        rectdf = rll[[1]]
        vccols = rll[[2]]
        ord = order(rectdf$x1)
        # Sort rectdf:
        rectdf = rectdf[ord,]
        vccols = vccols[ord]
        # Add grey - ubq rectangle:
        rectdf = rbind(c(x1=par()$usr[1], x2=rectdf$x1[1],
                         y1=par()$usr[3], y2=par()$usr[4]), rectdf)
        vccols = c('grey',vccols)
        rect(xleft=rectdf$x1, xright=rectdf$x2,
             ybottom=rectdf$y1, ytop=rectdf$y2,
             border=vccols, lwd=.75)
        par(xpd=NA)
        rect(xleft=rectdf$x1, xright=rectdf$x2, ybottom=par()$usr[3], 
             ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
             border='white', col=vccols, lwd=.25)
        par(xpd=FALSE)
    }
    par(mar=c(0.25, 0.25, 0, 0.25))
    clsord = dist.order[[set]][[1]]
    meta.image(enrichmat[clsord,1:5], colvals=colvals, cex=0, horiz=F, useRaster=TRUE)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    vbreaks = dist.order[[set]][[2]]
    abline(v=vbreaks,lty='dotted',lw=.5, col='darkgrey')
    # Plot counts:
    par(mar=c(0.25 + (!count.barplot) * 1.25, 0.25, 0, 0.25))
    if (count.barplot){
        plot.counts(counts, set, dist.order[[set]], cex.axis=0.75 - 0.2 * tiny)
    } else {
        # Break up the vccols:
        vcut = dist.order[[set]][[3]]
        rl = rle(vcut)$lengths
        vccols2 = c()
        for(i in 1:length(rl)){
            vccols2 = c(vccols2, rep(vccols[i], rl[i]))
        }
        plot.countprop(counts, set, dist.order[[set]], 
                       cex.axis=0.75 - 0.2 * tiny, cols=vccols2)
    }
    return(dist.order[[set]])
}




# ---------------------
# Frozen order:
# ---------------------
# NOTE: Not same cls as modules:
dford = read.delim('Annotation/bssid_order_frz20190326.tsv', header=F)
names(dford) <- c('id','name', 'cls')
dford = filter(dford, id %in% reducenames$id)
cellorder = as.character(dford$id)
rn = as.character(dford$id)
acut = dford$cls
dist.breaks = calc.breaks.acut(acut)
acut.nam = acut
names(acut.nam) = dford$id

# Labels in figure:
labels = meta[rn, 'GROUP']
faclabels = as.matrix(as.numeric(labels))
colset = as.character(rdcol$COLOR)
fix = rn[is.na(labels)]
print(length(fix))
# Label runs (that are not NONE):
lablist.orddist = label.runs(faclabels, labels, rdcol)

# --------------------------
# Ordered by centers matrix:
# --------------------------
method = 'ejaccard'
dt = dist(cmat, method)
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
# Labels in figure:
labelsmat = meta[rnmat, 'GROUP']
faclabels = as.matrix(as.numeric(labelsmat))
colset = as.character(rdcol$COLOR)
fix = rnmat[is.na(labelsmat)]
print(length(fix))
lablist.mat = label.runs(faclabels, labelsmat, rdcol)

# # ---------------------------------------
# # Ordered by groups, then centers matrix:
# # ---------------------------------------
# # Sort first by rdmp colors then by clusters
# ordgroup = order(labelsmat, decreasing=TRUE)
# rngroup = rnmat[ordgroup]
# labels = meta[rngroup, 'GROUP']
# faclabels = as.matrix(as.numeric(labels))
# colset = as.character(rdcol$COLOR)
# fix = rnmat[is.na(labels)]
# print(length(fix))
# acutgroup = as.numeric(labels)
# dist.breaks = calc.breaks.acut(acutgroup)
# acutgroup.nam = acutgroup
# names(acutgroup.nam) = rngroup
# lablist.group = label.runs(faclabels, labels, rdcol)

# ------------------------------------------
# V2: Ordered by groups, then centers matrix:
# ------------------------------------------
# Establish levels/ordering:
odf = rdcol  # Copy don't change rdcol
odf = rbind(odf, c('Multiple','grey80','Other'))
clvs = c('Tissue', 'Primary Cell', 'Immune', 'Neural', 'Stem-like', 'Other')
odf$category = factor(odf$category, levels=clvs)
# Sort alpha then group:
odf = odf[order(odf$GROUP),]
odf = odf[order(odf$category),]

# Sort first by odf/group ord then by clusters
labelsmat = as.character(labelsmat)
labelsmat = factor(labelsmat, levels=odf$GROUP)
ordgroup = order(labelsmat, decreasing=TRUE)
rngroup = rnmat[ordgroup]
labels = meta[rngroup, 'GROUP']
faclabels = as.matrix(as.numeric(labels))
colset = as.character(rdcol$COLOR)
fix = rnmat[is.na(labels)]
print(length(fix))
acutgroup = as.numeric(labels)
dist.breaks = calc.breaks.acut(acutgroup)
acutgroup.nam = acutgroup
names(acutgroup.nam) = rngroup
lablist.group = label.runs(faclabels, labels, rdcol)

