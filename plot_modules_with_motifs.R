#!/usr/bin/R
# -----------------------------------------------------
# Plot modules from CHMM/other clustering
# With MOTIFs 
# 1. Plot the centers matrix with motifs
# 2. Plot stats on motifs - amt of enrichment/depletion
# 3. Associate motifs with cell types by enrichment
# 4. Plot and order/cluster the adjacency matrix
#       for the epigenomes by motif co-enrichment
# -----------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'great_enrichment/auxfunctions_GREAT.R'))
library(png)

# Defaults:
motiffile='merged_enrichments_0_-99.tsv.gz'
filepref = 'cls_merge2_wH3K27ac100_300'
tagline = 'ChromHMM Enhancers'
imgdir = paste0(img, "clusters/") 

# Arguments:
args=(commandArgs(TRUE))
if (length(args) == 0) {
    stop("No arguments supplied: Need clusters filename")
} else {        
    motiffile = args[1]
    filepref = args[2]
    tagline = args[3]
    if (length(args) > 3){ 
        imgdir = args[4] 
    }
}

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){
    c(filepref, tagline, imgdir) }
source(paste0(bindir, 'load_modules_data.R'))

# Prefix:
cmd = paste('mkdir -p', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'clusters_motif_')
mpngdir = 'motif_data/logos2/'

# ---------------------------------
# Add motif enrichment information:
# ---------------------------------
# in our runs of test motif enrichment?
# ---------------------------------
# Load in enrichment files:
enrdf = read.delim(gzfile(motiffile), header=T)
print(length(unique(enrdf$Prefix)))
enrdf = enrdf[enrdf$compare == '+_.', c('Motif', 'Prefix', 'log2_enrich', 'log10_pval')]
wide = spread(enrdf[,c('Motif', 'Prefix', 'log2_enrich')], Motif, log2_enrich)
mmat = as.matrix(wide[,-1])
pwide = spread(enrdf[,c('Motif', 'Prefix', 'log10_pval')], Motif, log10_pval)
pmat = as.matrix(pwide[,-1])

# Strip the common prefix:
prefs = as.character(wide$Prefix)
preflen = min(sapply(prefs, nchar)) - 1 
numprefs = as.numeric(sapply(prefs, function(x){ substr(x, 61, 70) }))
rownames(mmat) = paste0('c',numprefs)
rownames(pmat) = paste0('c',numprefs)
# Keep only values crossing pvalue threshold
# mmat = mmat * (pmat >= 1.3)
print(dim(mmat))

# Motif information:
motfam = read.delim('Annotation/motifs-clust-names.txt', header=F, stringsAsFactors=F)
motmap = c()
for (i in 1:nrow(motfam)){
    members = strsplit(motfam[i,1], ";")[[1]]
    motmap = rbind(motmap, data.frame(motif = members, num=i-1,
                                      main = motfam[i,2], fam = motfam[i,3]))
}
rownames(motmap) = motmap$motif

# ---------------------------------
# Select only one motif per family:
# ---------------------------------
pick = 'max'
# pick = 'sum'
locmax = apply(abs(mmat), 2, max)
locsum = apply(abs(mmat), 2, sum)
if (pick == 'max'){
    maxdf = data.frame(topval=locmax, motif=names(locmax))
} else {
    maxdf = data.frame(topval=locsum, motif=names(locsum))
}
maxdf = merge(motmap, maxdf)
# keep.motifs = unlist(sapply(unique(motmap$num), function(i){
#                          subdf = maxdf[maxdf$num == i,]
#                          c(as.character(subdf[which.max(subdf$topval),'motif']),i) }))
# kmdf = ldply(unique(motmap$num), function(i){
kmdf = c()
for (i in unique(motmap$num)){
    subdf = maxdf[maxdf$num == i,]
    topm = as.character(subdf$motif[which.max(subdf$topval)])
    tval = subdf$topval[which.max(subdf$topval)]
    if (length(topm) > 0){ kmdf = rbind(kmdf, data.frame(motif=topm, num=i, topval=tval))}
}
# Keep motifs with at least logFC of cutoff (1.5)
cutoff = 1.5
keep.motifs = as.character(kmdf$motif[kmdf$topval >=cutoff])
redenr = mmat[,keep.motifs]
predenr = pmat[,keep.motifs]
print(dim(redenr))

# Numbers for manuscript:
# sum(apply(abs(mmat), 1, max) > 1.5) # 240 if all kept
sum(apply(abs(redenr), 1, max) > 1.5) # 202 of representative

# Motif order:
set = tagline
motifs = list()
keep.cls = list()
motifs[[set]] = redenr
keep.cls[[set]] = rownames(motifs[[set]])

# Variables for plot order:
rnorder = rngroup
acutnamed = acutgroup.nam
set = tagline 
lablist = lablist.group

# Diagonalize the motifs matrix:
dist.order = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, calc.only=TRUE)
clsord = dist.order[[1]]
vcut = dist.order[[3]]
mat = motifs[[set]]
ord = clsord[clsord %in% rownames(mat)]
ll = diag.mat2(abs(mat[ord,]))
col.motiford = ll[[2]]

plot.motifs <-  function(fmat, set, ordll, zlim=2, ablwd=1,
                       motord=NULL, labeled=TRUE, calc.only=FALSE){
    # GWAS PLOT HERE:
    ord = ordll[[1]]
    vbreaks = ordll[[2]]
    clscuts = ordll[[3]]
    subord = ord[ord %in% rownames(fmat)]
    gmat = fmat[subord,]
    if (!is.null(motord)){ gmat = gmat[,motord] }
    # Diagonalize, absolute val log2fc:
    tmat = abs(gmat)
    ll = diag.mat2(tmat, ratio=0.5)
    tmp = ll[[1]]
    cto = ll[[3]] # For breaks
    diagord = colnames(tmp)
    gmat = gmat[,diagord]
    vcut =  c(cto[cto ==0], clscuts[cto])
    hbreaks = calc.breaks.acut(vcut)
    gmat[gmat < -zlim] <- -zlim
    gmat[gmat > zlim] <- zlim
    # Threshold for plotting:
    if (!calc.only){
        image(gmat, axes=FALSE,col=rev(colrb), zlim=c(-zlim, zlim), useRaster=T)
        abline(h=par()$usr[3:4],lty=1,lw=0.5)
        abline(v=par()$usr[1:2],lty=1,lw=0.5)
        abline(v=vbreaks,lty='dotted',lw=ablwd, col='darkgrey')
        abline(h=hbreaks,lty='dotted',lw=ablwd, col='darkgrey')
    }
    return(list(diagord, vcut))
}

# Get order for clusters corresponding to motifs
ll = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed,
                  title=FALSE, subset=TRUE, calc.only=TRUE)
subcn = ll[[1]]
vbreaks = ll[[2]]
gll = plot.motifs(mmat, set, ll, motord=col.motiford, calc.only=TRUE)
motnam = gll[[1]]

# Plot:
png(paste0(imgpref, 'groupdist.png'), res=450,units='in',width=9,height=13)
ratio = nrow(redenr) / nrow(cmat) * 2
layout(matrix(c(1:12),4,3), heights=c(.6,8,8 * ratio, 1), widths=c(1.2,8,.5), TRUE)
par(yaxs="i")
par(xaxs="i")
dist.order = list()
# Metadata matrix:
par(mar=c(0.25, 6, 0, 0))
image(t(as.matrix(metamat) == ""), col='white', axes=F, useRaster=T)
metaclass = sapply(rev(colnames(metamat)), capitalize)
text(x=seq(0,1, length.out=ncol(metamat)),
     y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), 
     labels=metaclass, srt=90, adj=0, xpd=TRUE, cex=.7)
par(mar=c(.25, 6, 0, 0))
meta.image(metamat[rnorder,5:1], colvals=colvals, cex=0, horiz=T, useRaster=T)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
text(y=lablist[[1]],
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lablist[[2]], srt=0, adj=1, xpd=TRUE,cex=.8, col=lablist[[3]])
abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
# GWAS labels:
par(mar=c(0.25, 6, 0, 0))
tmpmat = motifs[[set]]
image(tmpmat, col='white', axes=F, useRaster=T)
text(y=seq(0, 1, length.out=ncol(tmpmat)),
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=motnam, srt=0, adj=1, xpd=TRUE,cex=.6)
# Add labels
plot.new()
par(mar=c(0.25, 0.25, 2, 0.25))
meta.image(enrichmat[subcn,5:1], colvals=colvals, cex=0, horiz=F, useRaster=T)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
abline(v=vbreaks,lty='dotted',lw=.5, col='darkgrey')
mtext(set, side=3, cex=1.3)
# Plot clusters and counts
set = tagline
par(mar=c(0.25, 0.25, 0, 0.25))
dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, 
                                 cls=acutnamed, title=FALSE, subset=TRUE)
rll = calc.breaks.rect(hcls=acutnamed, vcls=dist.order[[set]][[3]], colset)
clsrectdf = rll[[1]]
clsvccols = rll[[2]]
# Add grey - ubq rectangle:
clsrectdf = rbind(c(x1=par()$usr[1], x2=clsrectdf$x1[1],
                 y1=par()$usr[3], y2=par()$usr[4]), clsrectdf)
clsvccols = c('grey',clsvccols)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2,
     ybottom=clsrectdf$y1, ytop=clsrectdf$y2,
     border=clsvccols, lwd=1)
par(xpd=NA)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[3], 
     ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
par(xpd=FALSE)
abline(h=dist.breaks,lty='dotted', lw=.5, col='darkgrey')
# Plot motifs:
par(mar=c(0.25, 0.25, 0, 0.25))
gll = plot.motifs(mmat, set, dist.order[[set]], motord=col.motiford)
motnam = gll[[1]]
gcut = gll[[2]]
rll = calc.breaks.rect(hcls=gcut, vcls=dist.order[[set]][[3]], colset)
rectdf = rll[[1]]
vccols = rll[[2]]
# Add grey - ubq rectangle:
rectdf = rbind(c(x1=par()$usr[1], x2=rectdf$x1[1],
                 y1=par()$usr[3], y2=par()$usr[4]), rectdf)
vccols = c('grey',vccols)
rect(xleft=rectdf$x1, xright=rectdf$x2,
     ybottom=rectdf$y1, ytop=rectdf$y2,
     border=vccols, lwd=.75)
par(xpd=NA)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[3], 
     ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
par(xpd=FALSE)
par(mar=c(0.25, 0.25, 0, 0.25))
plot.counts(counts, set, dist.order[[set]])
plot.new()
# Availability:
par(mar=c(.25, 0, 0, 0.25))
avail = as.matrix(wm[main.marks, rnorder])
image(avail, axes=F, col=c('white', 'darkgrey'), useRaster=T)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
par(mar=c(.25, 0, 0, 0.25))
image(avail, axes=F, col='white', useRaster=T)
text(x=seq(0,1, length.out=length(main.marks)),
     y=par()$usr[4], 
     labels=main.marks, srt=90, adj=1, xpd=TRUE,cex=.5)
dev.off()



# ---------------------
# Plot only the motifs:
# ---------------------
# Pre-calculate:
ll = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed,
                  title=FALSE, subset=TRUE, calc.only=TRUE)
subcn = ll[[1]]
vbreaks = ll[[2]]
gll = plot.motifs(mmat, set, ll, motord=col.motiford, calc.only=TRUE)
motnam = gll[[1]]
gcut = gll[[2]]
gcol = c('grey', colset)[gcut + 1]

# Plot:
zlim = 1.5
# png(paste0(imgpref, 'groupdist_onlymotif_z', sub("\\.","_", zlim), '.png'), res=450, units='in', width=7.75,height=3)
pdf(paste0(imgpref, 'groupdist_onlymotif_z', sub("\\.","_", zlim), '.pdf'), width=7.75,height=3)
layout(matrix(c(1:2),1,2), widths=c(2, 8), TRUE)
par(yaxs="i")
par(xaxs="i")
dist.order = list()
sp = 0.15
# Motif labels:
par(mar=c(sp, 5, sp, sp))
tmpmat = motifs[[set]][, motnam]
image(tmpmat, col='white', axes=F, useRaster=T)
# Two columns of names + motif logos:
nb = 2
for (i in c(1:nb - 1)){
    # Select labels:
    fm = 1:length(motnam)
    fm = floor(fm / 8)
    fm = which((fm + 1) %% nb == i )
    # Adjust labels:
    box.pad = 0.017
    xold = seq(0, 1, length.out=ncol(tmpmat))[fm]
    xx = space.1d(xold, box.pad=box.pad)
    # xx = space.1d(xx, box.pad=box.pad)
    # Text:
    rt = i * 1.25 + 1
    text(y=xx, x=par()$usr[2]- rt*(par()$usr[2]-par()$usr[1]), 
         labels=motnam[fm], srt=0, adj=1, xpd=TRUE,cex=.35, col=gcol[fm])
    # Lines:
    xw = seq(0, 1, length.out=ncol(tmpmat))[fm]
    par(xpd=NA)
    if (i == 0){
        rx=c(0.005, 0.1, .4, .49, .5)
        x = par()$usr[2]-rx*(diff(par()$usr[1:2]))
        segments(x0=x[1], y0=xw, x1=x[2], y1=xw, col=gcol[fm], lwd=.5)
        segments(x0=x[2], y0=xw, x1=x[3], y1=xx, col=gcol[fm], lwd=.5)
        segments(x0=x[3], y0=xx, x1=x[4], y1=xx, col=gcol[fm], lwd=.5)
    } else {
        # Calculate centers of groups:
        sets = apply(matrix(c(rle(diff(fm))$lengths, 1), ncol=2, byrow=2), 1, sum)
        cs = c(0, cumsum(sets))
        f2 = rep(0, length(fm))  # assignments
        mx = rep(0, length(fm))  # centers
        nset = (length(cs)-1)
        for (j in 1:nset){
            f2[(cs[j]+1):cs[j+1]] = j
            sm = mean(xx[(cs[j]+1):cs[j+1]])
            # Adjust around the center:
            if (j == 1){ adj = -0.025 } else if (j == nset){ adj = 0.025 } else { adj = 0 }
            # mx[(cs[j]+1):cs[j+1]] = sm + adj
            rnj = (cs[j]+1):cs[j+1]
            rnj = rnj - min(rnj)
            pad = 0.018
            mx[(cs[j]+1):cs[j+1]] = sm + adj + (pad * (rnj - 1)/(max(rnj) - 1) - pad/2)
            # Adjust to space:
        }
        rx=c(0.005, 0.1, .4, 1.5, 1.65, 1.74, 1.75)
        x = par()$usr[2]-rx*(diff(par()$usr[1:2]))
        segments(x0=x[1], y0=xw, x1=x[2], y1=xw, col=gcol[fm], lwd=.5)
        segments(x0=x[2], y0=xw, x1=x[3], y1=mx, col=gcol[fm], lwd=.5)
        segments(x0=x[3], y0=mx, x1=x[4], y1=mx, col=gcol[fm], lwd=.5)
        segments(x0=x[4], y0=mx, x1=x[5], y1=xx, col=gcol[fm], lwd=.5)
        segments(x0=x[5], y0=xx, x1=x[6], y1=xx, col=gcol[fm], lwd=.5)
    }
    # Add the logos:
    for (fj in 1:length(fm)){
        j = fm[fj] # Motif ind
        yj = xx[fj] # Y-coord
        mj = motnam[j]
        img <- try(readPNG(paste0(mpngdir, mj, '.png')))
        xj = par()$usr[2]- (rt - .5)*(par()$usr[2]-par()$usr[1])
        ysize = .008
        xsize = .5
        # Image is 600 x 800
        if (class(img) != 'try-error'){
            rasterImage(img, xleft=xj - xsize, xright=xj,
                        ybottom=yj-ysize, ytop=yj + ysize)
            rect(xleft=xj, xright=xj + ysize * 2, 
                 ybottom=yj-ysize, ytop=yj + ysize, 
                 col=gcol[j], border=NA, lwd=0.5)
        }
    }
    par(xpd=FALSE)
}
# Calculate clusters and counts
set = tagline
dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, 
                                 cls=acutnamed, title=FALSE, subset=TRUE, calc.only=TRUE)
rll = calc.breaks.rect(hcls=acutnamed, vcls=dist.order[[set]][[3]], colset)
clsrectdf = rll[[1]]
clsvccols = rll[[2]]
clsrectdf = rbind(c(x1=par()$usr[1], x2=clsrectdf$x1[1],
                    y1=par()$usr[3], y2=par()$usr[4]), clsrectdf)
clsvccols = c('grey',clsvccols)
# Plot motifs:
par(mar=c(sp, 0, sp, sp))
gll = plot.motifs(mmat, set, dist.order[[set]], motord=col.motiford, ablwd=.5, zlim=zlim)
motnam = gll[[1]]
gcut = gll[[2]]
rll = calc.breaks.rect(hcls=gcut, vcls=dist.order[[set]][[3]], colset)
rectdf = rll[[1]]
vccols = rll[[2]]
# Add grey - ubq rectangle:
rectdf = rbind(c(x1=par()$usr[1], x2=rectdf$x1[1],
                 y1=par()$usr[3], y2=par()$usr[4]), rectdf)
vccols = c('grey',vccols)
rect(xleft=rectdf$x1, xright=rectdf$x2,
     ybottom=rectdf$y1, ytop=rectdf$y2,
     border=vccols, lwd=.75)
par(xpd=NA)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[3], 
     ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[4], 
     ytop=par()$usr[4]+0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
par(xpd=FALSE)
dev.off()




add.motifs = function(fm, motnam, flip=FALSE, motcex=0.35){
    box.pad = 0.02
    xold = seq(0, 1, length.out=ncol(tmpmat))[fm]
    xx = space.1d(xold, box.pad=box.pad)
    xx = space.1d(xx, box.pad=box.pad)
    # Text:
    rt = 1
    if (flip){
        xat = par()$usr[1]+ rt*(par()$usr[2]-par()$usr[1])
    } else { 
        xat = par()$usr[2]- rt*(par()$usr[2]-par()$usr[1])
    }
    text(y=xx, x=xat, labels=motnam[fm], 
         srt=0, adj=1 - flip, xpd=TRUE,cex=motcex, col=gcol[fm])
    # Lines:
    xw = seq(0, 1, length.out=ncol(tmpmat))[fm]
    par(xpd=NA)
    rx=c(0.005, 0.1, .4, .49, .5)
    if(flip){
        x = par()$usr[1]+rx*(diff(par()$usr[1:2]))
    } else {
        x = par()$usr[2]-rx*(diff(par()$usr[1:2]))
    }
    segments(x0=x[1], y0=xw, x1=x[2], y1=xw, col=gcol[fm], lwd=.5)
    segments(x0=x[2], y0=xw, x1=x[3], y1=xx, col=gcol[fm], lwd=.5)
    segments(x0=x[3], y0=xx, x1=x[4], y1=xx, col=gcol[fm], lwd=.5)
    # Add the logos:
    for (fj in 1:length(fm)){
        j = fm[fj] # Motif ind
        yj = xx[fj] # Y-coord
        mj = motnam[j]
        img <- try(readPNG(paste0(mpngdir, mj, '.png')))
        if (flip){
            xj = par()$usr[1] + (rt - .5)*(par()$usr[2]-par()$usr[1])
        } else {
            xj = par()$usr[2]- (rt - .5)*(par()$usr[2]-par()$usr[1])
        }
        ysize = .008 * (motcex / 0.35)  # Scale with text
        xsize = .5
        # Image is 600 x 800
        if (class(img) != 'try-error'){
            if (flip){
                xat = c(xj, xj + xsize)
            } else {
                xat = c(xj - xsize, xj)
            }
            rasterImage(img, xleft=xat[1], xright=xat[2],
                        ybottom=yj-ysize, ytop=yj + ysize)
            rect(xleft=xj, xright=xj + ysize * (2 - 4 * flip), 
                 ybottom=yj-ysize, ytop=yj + ysize, 
                 col=gcol[j], border=NA, lwd=0.5)
        }
    }
    par(xpd=FALSE)
}

# Plot:
zlim = 1.5
# png(paste0(imgpref, 'groupdist_onlymotif_z', sub("\\.","_", zlim), '.png'), res=450, units='in', width=7.75,height=3)
pdf(paste0(imgpref, 'groupdist_onlymotif_bothsides_z', sub("\\.","_", zlim), '.pdf'), width=8.75,height=3)
layout(matrix(c(1:3),1,3), widths=c(2, 8, 2), TRUE)
par(yaxs="i")
par(xaxs="i")
dist.order = list()
sp = 0.15
motcuts = rev(c(9, 5, 6, 6, 8, 6, 6, 6, 5, 7, 8, 8, 7))
motside = unlist(sapply(1:length(motcuts), function(i){ rep(i%%2, motcuts[i])}))
motside = motside[1:length(motnam)]
# Motif labels:
par(mar=c(sp, 5, sp, sp))
tmpmat = motifs[[set]][, motnam]
image(tmpmat, col='white', axes=F, useRaster=T)
# Single column motifs:
fm = which(motside == 1)
add.motifs(fm, motnam, flip=FALSE, motcex=.45)
# Calculate clusters and counts
set = tagline
dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, 
                                 cls=acutnamed, title=FALSE, subset=TRUE, calc.only=TRUE)
rll = calc.breaks.rect(hcls=acutnamed, vcls=dist.order[[set]][[3]], colset)
clsrectdf = rll[[1]]
clsvccols = rll[[2]]
clsrectdf = rbind(c(x1=par()$usr[1], x2=clsrectdf$x1[1],
                    y1=par()$usr[3], y2=par()$usr[4]), clsrectdf)
clsvccols = c('grey',clsvccols)
# Plot motifs:
par(mar=c(sp, 0, sp, 0))
gll = plot.motifs(mmat, set, dist.order[[set]], motord=col.motiford, ablwd=.5, zlim=zlim)
motnam = gll[[1]]
gcut = gll[[2]]
rll = calc.breaks.rect(hcls=gcut, vcls=dist.order[[set]][[3]], colset)
rectdf = rll[[1]]
vccols = rll[[2]]
# Add grey - ubq rectangle:
rectdf = rbind(c(x1=par()$usr[1], x2=rectdf$x1[1],
                 y1=par()$usr[3], y2=par()$usr[4]), rectdf)
vccols = c('grey',vccols)
rect(xleft=rectdf$x1, xright=rectdf$x2,
     ybottom=rectdf$y1, ytop=rectdf$y2,
     border=vccols, lwd=.75)
par(xpd=NA)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[2], 
     xleft=par()$usr[2]+0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[3], 
     ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[4], 
     ytop=par()$usr[4]+0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
par(xpd=FALSE)
par(mar=c(sp, sp, sp, 5))
tmpmat = motifs[[set]][, motnam]
image(tmpmat, col='white', axes=F, useRaster=T)
fm = which(motside == 0)
add.motifs(fm, motnam, flip=TRUE, motcex=.45)
dev.off()


# Statistics on the motifs (for main text):
mat = abs(redenr)
mmarg = apply(mat > 1.5, 2, sum)
# mmarg = apply(mat > 1, 2, sum)
print(sort(mmarg) / 300)
kind = names(which(mmarg /300 < 0.05))
mean(mmarg[kind])

cmarg = apply(mat > 1.5, 1, sum)
print(sum(cmarg > 0))
cmarg = apply(mat > 1, 1, sum)
print(sum(cmarg > 0))

n = motmap$num[motmap$motif == 'HNF1A_4']
motmap[motmap$num == n, ]
# RFX2_3: RFX1,2,3,4,5
# HNF1A4: HNF1;HNF1A1,2,3;HNF1B1,2,3



