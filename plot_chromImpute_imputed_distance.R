#!/usr/bin/R
# ---------------------------------------
# Plot the imputed, observed, and 
# mix of imputed/observed matrices
# Plot MDS
# Validate imputed diff
# Look for consistency in clades
# Plot dendrograms colored by tissue type
# Identify consistent breakpoints
# ---------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
library(ggplot2)
library(viridis)
library(ggrepel)
library(dendextend)
library(uwot)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
options(scipen=45) # So we dont get issues writing integers into bedfiles

# Load specific distance matrices:
fixedstate = FALSE
if (fixedstate){
    source(paste0(bindir, 'load_region_distance_matrices.R'))
    setprefix = paste0('region_',nregions,'_distances_')
} else { 
    source(paste0(bindir, 'load_distance_matrices.R'))
    setprefix = 'distances_'
}

today <- format(Sys.time(), "%m%d%y")
imgdir = paste0(img, "imp_distance/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, setprefix)

# ------------------------------
# Reorder fused distance matrix:
# ------------------------------
dt <- as.dist(full)
ht <- hclust(dt, method='ward.D')
cocl <- order.optimal(dt, ht$merge)$order
rn <- names(cocl)[cocl]
full = full[rn,rn]
# Generate breaks:
NCLUST = 20
breaks = calc.breaks(ht, NCLUST, cocl)
acut <- cutree(ht, NCLUST)[cocl]

# Save this order:
rdf = merge(data.frame(id=rn), map)
rdf$id = factor(rdf$id, levels=rn)
rdf = rdf[order(rdf$id),]
rdf$cls = acut
write.table(rdf, paste0('Annotation/',setprefix, 'bssid_order_20190219.tsv'), row.names=F, col.names=F, quote=F, sep="\t")

# Accessory
tmp <- matrix(NA, nrow=N, ncol=N, dimnames=list(rn,rn))
palette = viridis(100)
palette = colryb

# Make image:
# png(paste0(imgpref,'matrix_comparison_all.png'),res=300,units='in',width=9,height=4.5)
pdf(paste0(imgpref,'matrix_comparison_all.pdf'),width=9,height=4.5)
layout(matrix(c(1,1,2,3,1,1,4:(2 * (NMARKS + 1) + 1)), 4, 8),
       heights=rep(1,4), widths=rep(1, 8), TRUE)
sp=0.1
par(mar=c(sp,sp, 2, sp))
plot.cov(full, clamp=TRUE, palette=palette, breaks=breaks)
mtext('Fused', side=3, cex=1)
# Add other images:
for (mark in marks[-length(marks)]){ 
    # Observed:
    par(mar=c(sp,sp, 2, sp))
    mat = obsll[[mark]]
    subrn = rn[rn %in% colnames(mat)]
    mat2 = tmp
    mat2[subrn,subrn] = mat[subrn, subrn]
    plot.cov(mat2, clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey')
    mtext(mark, side=3, cex=1)
    # Imputed:
    par(mar=rep(sp,4))
    mat = ll[[mark]]
    # mat = clamp.mat(mat)
    subrn = rn[rn %in% colnames(mat)]
    mat2 = tmp
    mat2[subrn,subrn] = mat[subrn, subrn]
    plot.cov(mat2, clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey')
}
dev.off()

# ------------------------
# Plot matrix with colors:
# ------------------------
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
fix = rn[is.na(labels)]
print(length(fix))
# Label runs (that are not NONE):
lablist = label.runs(faclabels, labels, rdcol)

add.lablist.text = function(lablist, box.pad=0.02, rx=c(0.1,0.15, 0.3,0.35), small=TRUE, cex=.8){
    par(xpd=TRUE)
    xx = space.1d(lablist[[1]], box.pad=box.pad, lim=c(0 + 3e-3, 1 - 3e-3))
    xx = space.1d(xx, box.pad=box.pad, lim=c(0 + 3e-3, 1 - 3e-3))
    x = par()$usr[1]-rx*(diff(par()$usr[1:2]))
    text(y=xx, x=x[4], labels=lablist[[2]],
         srt=0, adj=1, xpd=TRUE, cex=cex, col=lablist[[3]])
    segments(x0=x[3], y0=xx, x1=x[2], y1=lablist[[1]], col=lablist[[3]])
    segments(x0=x[1], y0=lablist[[1]], x1=x[2], y1=lablist[[1]], col=lablist[[3]])
    par(xpd=FALSE)
}

# Matrix with roadmap colors
png(paste0(imgpref,'matrix_rdcol.png'),res=300,units='in',width=12,height=12)
layout(matrix(c(1,2,3,4),2,2), widths=c(1.5,8), heights=c(2,8), TRUE)
par(yaxs="i")
par(xaxs="i")
# Empty square:
par(mar=c(0.5, 0.5, 0.5, 0.5))
plot.new()
par(mar=c(0.5, 10, 0, 0.05))
image(t(faclabels), axes=F, col=colset, useRaster=TRUE)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
add.lablist.text(lablist, box.pad=0.01, rx=c(0.1,0.3, 0.7,0.8))
par(mar=c(0, 0.05, 2, 0.5))
plot(ht, main='', col='black',lwd=1,lty=1,hang=-1,labels=FALSE,axes=F,ylab='')
abline(v=breaks * nrow(full),lty=1,lwd=.5, col='darkblue')
par(mar=c(0.5, 0.05, 0, 0.5))
plot.cov(full, clamp=TRUE, palette=palette, breaks=breaks)
dev.off()

tmp = full
tmp = tmp * NA
# Matrix with roadmap colors + all others + availability + metadata
# NOTE: Same as above with avail.
small=TRUE
# png(paste0(imgpref,'matrix_all_rdcol_withavail.png'),res=300,units='in',width=8.5,height=4.5)
pdf(paste0(imgpref,'matrix_all_rdcol_withavail.pdf'),width=8.5,height=4.5)
layout(matrix(c(1,rep(2,4), 3,rep(4,4), 3, rep(4,4), 5, rep(6,4), 7:26),5,8), 
       widths=c(1.25 + 1.25 * small,4,4, .75 + 0.25 * small,2,2,2,2), heights=rep(2,5), TRUE)
sp = 0.15
par(yaxs="i")
par(xaxs="i")
par(mar=c(sp, 6, sp, sp))
image(t(as.matrix(metamat) == ""), col='white', axes=F, useRaster=TRUE)
text(x=seq(0,1, length.out=ncol(metamat)),
     y=par()$usr[3]+0.005*(par()$usr[4]-par()$usr[3]), 
     labels=sapply(rev(colnames(metamat)), capitalize),
     srt=90, adj=0, xpd=TRUE,cex=1 - .5 * small)
par(mar=c(sp, 6, 0, 0.15))
meta.image(metamat[rownames(full),5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
abline(h=breaks,lty='dotted',lwd=.25, col='black')
box(lwd=0.5)
add.lablist.text(lablist, box.pad=0.016, rx=c(0.05,0.1, 0.3,0.35), cex=.6)
par(mar=c(0, 0.05, sp, sp))
plot(ht, main='', col='black',lwd=.5,lty=1,hang=-1,labels=FALSE,axes=F,ylab='', xlab='')
abline(v=breaks * nrow(full),lty=1,lwd=.25, col='black')
par(mar=c(sp, 0.05, 0, sp))
plot.cov(full, clamp=TRUE, palette=palette, breaks=breaks, blty=1)
# Add other marks:
par(mar=c(sp, sp, sp, sp))
avail = as.matrix(wm[main.marks, rn])
image(avail, axes=F, col='white', useRaster=TRUE)
text(x=seq(0,1, length.out=length(main.marks)),
     y=par()$usr[3]+0.01*(par()$usr[4]-par()$usr[3]), 
     labels=main.marks, srt=90, adj=0, xpd=TRUE,cex=.6 - 0.2 * small)
par(mar=c(sp, sp, 0, sp))
image(avail, yaxt='n',xaxt='n', col=c('white','darkgrey'), useRaster=TRUE)
abline(h=breaks,lty='dotted',lwd=.25, col='black')
for (i in 1:(NMAIN + 1)){ 
    mark = c('DNase-seq', mainmarks)[i]
    if (i %% 2 == 1){ 
        plot.new()
        par(mar=c(sp, sp, 0, sp))
    } else {
        par(mar=c(sp, sp, 2 - 1 * small, sp))
    }
    # Observed:
    mat = obsll[[mark]]
    subrn = rn[rn %in% colnames(mat)]
    mat2 = tmp
    mat2[subrn,subrn] = mat[subrn, subrn]
    plot.cov(mat2, clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey')
    mtext(mark, side=3, cex=1.3 - .7 * small)
    # Imputed:
    par(mar=rep(sp,4))
    mat = ll[[mark]]
    subrn = rn[rn %in% colnames(mat)]
    mat2 = tmp
    mat2[subrn,subrn] = mat[subrn, subrn]
    plot.cov(mat2, clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey')
}
dev.off()


tmp = full
tmp = tmp * NA
# Matrix with roadmap colors + all others + availability + metadata
# NOTE: Same as above with avail.
small=TRUE
zlim = c(0, max(full, na.rm=T))
ozlim = c(0,.1)
for (i in 1:NMAIN){ 
    mark = mainmarks[i]
    zlim[2] = max(max(ll[[mark]], na.rm=T), zlim[2])
    print(zlim)
    ozlim[2] = max(max(obsll[[mark]], na.rm=T), ozlim[2])
    print(ozlim)
}
zlim =ozlim
palette = rev(viridis(100))
png(paste0(imgpref,'matrix_all_rdcol_withavail_samescale.png'),res=300,units='in',width=8.25,height=4.75)
# pdf(paste0(imgpref,'matrix_all_rdcol_withavail_samescale.pdf'),width=8.25,height=4.75)
layout(matrix(c(1,rep(2,4), 3,rep(4,4), 3, rep(4,4), 5, rep(6,4), 7:21),5,7), 
       widths=c(1.25 + 1.25 * small,4,4, .75 + 0.25 * small,2,2,2), heights=rep(2,5), TRUE)
sp = 0.15
par(yaxs="i")
par(xaxs="i")
par(mar=c(sp, 6, sp, sp))
image(t(as.matrix(metamat) == ""), col='white', axes=F, useRaster=TRUE)
text(x=seq(0,1, length.out=ncol(metamat)),
     y=par()$usr[3]+0.005*(par()$usr[4]-par()$usr[3]), 
     labels=sapply(rev(colnames(metamat)), capitalize),
     srt=90, adj=0, xpd=TRUE,cex=1 - .5 * small)
par(mar=c(sp, 6, 0, 0.15))
meta.image(metamat[rownames(full),5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
add.lablist.text(lablist, box.pad=0.016, rx=c(0.05,0.1, 0.3,0.35), cex=.6)
par(mar=c(0, 0.05, sp, sp))
plot(ht, main='', col='black',lwd=1,lty=1,hang=-1,labels=FALSE,axes=F,ylab='', xlab='')
abline(v=breaks * nrow(full),lty=1,lwd=.5, col='black')
par(mar=c(sp, 0.05, 0, sp))
plot.cov(full, clamp=TRUE, palette=palette, breaks=breaks, zlim=zlim)
# Add other marks:
par(mar=c(sp, sp, sp, sp))
avail = as.matrix(wm[main.marks, rn])
image(avail, axes=F, col='white', useRaster=TRUE)
text(x=seq(0,1, length.out=length(main.marks)),
     y=par()$usr[3]+0.01*(par()$usr[4]-par()$usr[3]), 
     labels=main.marks, srt=90, adj=0, xpd=TRUE,cex=.6 - 0.2 * small)
par(mar=c(sp, sp, 0, sp))
image(avail, yaxt='n',xaxt='n', col=c('white','darkgrey'), useRaster=TRUE)
abline(h=breaks,lty=1,lwd=.5, col='black')
for (i in 1:NMAIN){ 
    mark = mainmarks[i]
    if (i %% 2 == 1){ 
        plot.new()
        par(mar=c(sp, sp, 0, sp))
    } else {
        par(mar=c(sp, sp, 2 - 1 * small, sp))
    }
    # Observed:
    mat = obsll[[mark]]
    subrn = rn[rn %in% colnames(mat)]
    mat2 = tmp
    mat2[subrn,subrn] = mat[subrn, subrn]
    plot.cov(mat2, clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey', zlim=ozlim)
    mtext(mark, side=3, cex=1.3 - .7 * small)
    # Imputed:
    par(mar=rep(sp,4))
    mat = ll[[mark]]
    subrn = rn[rn %in% colnames(mat)]
    mat2 = tmp
    mat2[subrn,subrn] = mat[subrn, subrn]
    plot.cov(mat2, clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey', zlim=zlim)
}
dev.off()


# ----------------------------
# Update for all imputed marks
# ----------------------------
tmp = full
tmp = tmp * NA
# Matrix with roadmap colors + all others + availability + metadata
# NOTE: Same as above with avail.
# png(paste0(imgpref,'matrix_all_rdcol_withavail_allmarks.png'),res=300,units='in',width=25,height=11)
# png(paste0(imgpref,'matrix_all_rdcol_withavail_allmarks.png'),res=300,units='in',width=8.25,height=3.5)
pdf(paste0(imgpref,'matrix_all_rdcol_withavail_allmarks.pdf'),width=8.25,height=3.5)
layout(matrix(c(1,rep(2,4), 3,rep(4,4), 3, rep(4,4), 5, rep(6,4), 7:41),5,11), 
       widths=c(1.75,4,4, .75,rep(2, 8)), heights=rep(2,5), TRUE)
par(yaxs="i")
par(xaxs="i")
sp = 0.15
par(mar=c(sp, 3, sp,sp))
image(t(as.matrix(metamat) == ""), col='white', axes=F, useRaster=TRUE)
text(x=seq(0,1, length.out=ncol(metamat)),
     y=par()$usr[3]+0.01*(par()$usr[4]-par()$usr[3]), 
     labels=rev(colnames(metamat)), srt=90, adj=0, xpd=TRUE,cex=.3)
par(mar=c(sp, 3, 0, sp))
meta.image(metamat[rownames(full),5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
text(y=lablist[[1]],
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lablist[[2]], srt=0, adj=1, xpd=TRUE,cex=.4, col=lablist[[3]])
par(mar=c(0, sp, 2, sp))
plot(ht, main='', col='black',lwd=1,lty=1,hang=-1,labels=FALSE,axes=F,ylab='')
mtext('Fused', side=3, cex=.5)
abline(v=breaks * nrow(full),lty=1,lwd=.5, col='darkblue')
par(mar=c(sp, sp, 0, sp))
plot.cov(full, clamp=TRUE, palette=palette, breaks=breaks)
# Add other marks:
par(mar=c(0, sp, 2, sp))
avail = as.matrix(wm[main.marks, rn])
image(avail, axes=F, col='white', useRaster=TRUE)
text(x=seq(0,1, length.out=length(main.marks)),
     y=par()$usr[3]+0.01*(par()$usr[4]-par()$usr[3]), 
     labels=main.marks, srt=90, adj=0, xpd=TRUE,cex=.3)
par(mar=c(sp, sp, 0, sp))
image(avail, yaxt='n',xaxt='n', col=c('white','darkgrey'), useRaster=TRUE)
abline(h=breaks,lty=1,lwd=.5, col='black')
for (i in 1:NMARKS){ 
    mark = marks[i]
    if (i %% 2 == 1){ 
        plot.new()
        par(mar=c(sp, sp, 0, sp))
    } else {
        par(mar=c(sp, sp, .5, sp))
    }
    # Observed:
    mat = obsll[[mark]]
    subrn = rn[rn %in% colnames(mat)]
    mat2 = tmp
    mat2[subrn,subrn] = mat[subrn, subrn]
    plot.cov(mat2, clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey')
    mtext(mark, side=3, cex=.4)
    # Imputed:
    par(mar=rep(sp,4))
    mat = ll[[mark]]
    subrn = rn[rn %in% colnames(mat)]
    mat2 = tmp
    mat2[subrn,subrn] = mat[subrn, subrn]
    plot.cov(mat2, clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey')
}
dev.off()

# ------------------------------------------
# Sort first by rdmp colors then by clusters
# ------------------------------------------
rn = rn[order(labels, decreasing=TRUE)]
full = full[rn,rn]
labels = meta[rn, 'GROUP']
faclabels = as.matrix(as.numeric(labels))
colset = as.character(rdcol$COLOR)
# Loc labels. Adjust vals if very close to each other.
ldf = data.frame(label=labels, quant=1:length(labels) * 1.0/ length(labels) )
labeldf = aggregate(quant ~ label, ldf, mean)

# Matrix with roadmap colors + all others:
png(paste0(imgpref,'matrix_all_BYrdcol.png'),res=300,units='in',width=18,height=12)
layout(matrix(c(rep(1,4),rep(2,4), rep(2,4), 3:14),4,6), widths=c(1,4,4, 2,2,2), 
       heights=c(2,2,2,2), TRUE)
par(yaxs="i")
par(xaxs="i")
par(mar=c(0.5, 6, 3, 0.05))
meta.image(metamat[rownames(full),5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
text(x=seq(0,1, length.out=ncol(metamat)),
     y=par()$usr[4]+0.0025*(par()$usr[4]-par()$usr[3]), 
     labels=colnames(metamat), srt=90, adj=0, xpd=TRUE,cex=.9)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
text(y=labeldf$quant,
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=labeldf$label, srt=0, adj=1, xpd=TRUE,cex=.8)
par(mar=c(0.5, 0.05, 3, 0.5))
plot.cov(full, clamp=TRUE, palette=palette, breaks=breaks)
mtext('Fused', side=3, cex=1.4, line=.5)
# Add other marks:
for (i in 1:NMAIN){ 
    mark = mainmarks[i]
    par(mar=c(0.5, 0.5, 3, 0.5))
    # Observed:
    mat = obsll[[mark]]
    subrn = rn[rn %in% colnames(mat)]
    mat2 = tmp
    mat2[subrn,subrn] = mat[subrn, subrn]
    plot.cov(mat2, clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey')
    mtext(mark, side=3, cex=1.4, line=.5)
    # Imputed:
    par(mar=rep(0.5,4))
    mat = ll[[mark]]
    subrn = rn[rn %in% colnames(mat)]
    mat2 = tmp
    mat2[subrn,subrn] = mat[subrn, subrn]
    plot.cov(mat2, clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey')
}
dev.off()

