#!/usr/bin/R
# -----------------------------------------------------
# Plot modules from CHMM/other clustering with Motifs 
# 1. Plot the centers matrix with gwas
# 2. Plot stats on gwas - amt of enrichment/depletion
# 3. Associate gwas with cell types by enrichment
# 4. Plot and order/cluster the adjacency matrix
#       for the epigenomes by gwas co-enrichment
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
library(viridis)
library(tidyr)
library(dplyr)

# Defaults:
motiffile='merged_enrichments_0_-99.tsv'
filepref = 'cls_merge2_wH3K27ac100_raw'
tagline = 'ChromHMM Enhancers (on Epigenomes)'
imgdir = paste0(img, "clusters/") 

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need clusters filename")
} else {        
    motiffile = args[1]
    filepref = args[2]
    tagline = args[3]
    if (length(args) > 3){ 
        imgdir = args[4] 
    }
}

# Prefix / plot directory:
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'motifs_',filepref,'_')
print(paste("[STATUS] Plotting under prefix", imgpref))


# -----------------
# Load motifs data:
# -----------------
# Load in enrichment files:
enrdf = read.delim(motiffile, header=T)
print(length(unique(enrdf$Prefix)))
enrdf = enrdf[enrdf$compare == '+_.', c('Motif', 'Prefix', 'log2_enrich')]
wide = spread(enrdf, Motif, log2_enrich)
mmat = as.matrix(wide[,-1])

# Strip the common prefix:
prefs = as.character(wide$Prefix)
preflen = min(sapply(prefs, nchar)) - 1 
numprefs = as.numeric(sapply(prefs, function(x){ substr(x, 61, 70) }))
rownames(mmat) = paste0('c',numprefs)
print(dim(mmat))

# Motif information:
motfam = read.delim('Annotation/motifs-clust-names.txt', header=F, stringsAsFactors=F)
motmap = c()
for (i in 1:nrow(motfam)){
    members = strsplit(motfam[i,1], ";")[[1]]
    motmap = rbind(motmap, data.frame(motif = members, num=i,
                                      main = motfam[i,2], fam = motfam[i,3]))
}
rownames(motmap) = motmap$motif

# Select only one motif per family:
#	sum(colnames(mmat) %in% motmap$motif) / ncol(mmat) * 100
locmax = apply(abs(mmat),2, max)
maxdf = data.frame(topval=locmax, motif=names(locmax))
maxdf = merge(motmap, maxdf)
keep.motifs = unlist(sapply(unique(motmap$num), function(i){
                         subdf = maxdf[maxdf$num == i,]
                         as.character(subdf[which.max(subdf$topval),'motif']) }))
keep.motifs = keep.motifs[locmax[keep.motifs] >= 1.5]
redenr = mmat[,keep.motifs]
print(dim(redenr))

redenr[redenr < -2] = 2
redenr[redenr > 2] = 2
image(redenr, col=rev(colrb), zlim=c(-2,2))


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
main = c()
mat = motifs[[set]]
ord = clsord[clsord %in% rownames(mat)]
col.motiford = diag.mat(mat[ord,])[[2]]

# Ratios:
widths = unlist(lapply(motifs, nrow)) / nrow(redenr) * 10
clsord = dist.order[[1]]
cn = clsord[clsord %in% rownames(mmat)]
range(mmat)
image(t(mmat[cn,]), col=rev(colrb), zlim=c(-2,2))


# Get order for clusters corresponding to motifs
ll = plot.centers(centlist, set, rnorder,
                     counts=counts, cls=acutnamed,
                     title=FALSE, subset=TRUE,
                     calc.only=TRUE)
subcn = ll[[1]]
vbreaks = ll[[2]]



























# --------------------------------
# Add gwas enrichment information:
# --------------------------------
gwdf = read.delim(gwasfile, header=F)
names(gwdf) = c('pvalue','cluster','pmid','trait',
              'counthit','countall','fold')
gwdf$set = tagline
gwdf$pmt = paste0(gwdf$pmid, '_', gwdf$trait)
gwdf$cls = paste0('c', gwdf$cluster)
gwdf$logpval = -log10(gwdf$pvalue)
gwlong = aggregate(logpval ~ cls + pmt, gwdf, max)
wide = spread(gwlong, pmt, logpval, fill=0)
gwmat = as.matrix(wide[,-1])
rownames(gwmat) = wide$cls
gwmat[gwmat < 1] = 0

# Choose # gwas to show:
SHOWGWAS=350
# gwasmarg = sort(apply(gwmat, 2, sum), decreasing=T)
gwasmarg = sort(apply(gwmat, 2, max), decreasing=T)
keep.studies = names(head(gwasmarg, SHOWGWAS))
zmax = 5
zmin=2

# Order the top studies:
r2 = reord(t(gwmat[, keep.studies]) > zmin, 'Jaccard')
studyord = rownames(r2)
r3 = reord(gwmat[, keep.studies] > zmin, 'eJaccard')
roword = rownames(r3)
finalmat = gwmat[roword, studyord]

# Threshold for plotting:
gmat = finalmat
gmat[gmat > zmax] <- zmax
gmat[gmat < zmin] <- 0

# Diagonalize, ignore all <2
tmat = gwmat
tmat[tmat < zmin] = 0
ll = diag.mat(tmat[roword,keep.studies])
tmp = ll[[1]]
cto = ll[[3]] # For breaks
tmp[tmp > zmax] <- zmax

# Plot gwas alone 
png(paste0(imgpref,'top',SHOWGWAS,'_alone.png'),res=450,units='in',width=17,height=15)
par(mar=c(0.25, 15, .25, .25))
image(gmat, axes=FALSE,col=colred, zlim=c(zmin, zmax))
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
# abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
rows = (seq(0, ncol(gmat), 10) - .5) / (ncol(gmat) - 1)
cols = (seq(0, nrow(gmat), 20) - .5) / (nrow(gmat) - 1)
abline(h=rows,lty='dotted',lw=1, col='darkgrey')
abline(v=cols,lty='dotted',lw=1, col='darkgrey')
text(y=seq(0,1, length.out=length(studyord)),
     x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=studyord, srt=0, adj=1, xpd=TRUE,cex=.35)
dev.off()


png(paste0(imgpref,'top',SHOWGWAS,'_alone_diag.png'),res=450,units='in',width=17,height=15)
par(mar=c(0.25, 15, .25, .25))
image(tmp, axes=FALSE,col=colred, zlim=c(zmin, zmax))
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
# abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
rows = (seq(0, ncol(gmat), 10) - .5) / (ncol(gmat) - 1)
cols = (seq(0, nrow(gmat), 20) - .5) / (nrow(gmat) - 1)
abline(h=rows,lty='dotted',lw=1, col='darkgrey')
abline(v=cols,lty='dotted',lw=1, col='darkgrey')
text(y=seq(0,1, length.out=ncol(tmp)),
     x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=colnames(tmp), srt=0, adj=1, xpd=TRUE,cex=.35)
dev.off()


# ------------------------------
# After plot alone, try matched:
# ------------------------------
namesfile = paste0(filepref, '_names.tsv')
epinames = scan(namesfile,'c')
is.epi = length(grep("BSS",epinames)) == nrow(gwmat)
if (is.epi){
    print("Is enrichment on epigenomes, plotting with epigenomes")
    clsn = paste0('c',1:length(epinames) - 1)
    names(epinames) = clsn

    # Map:
    epimat = gwmat
    rownames(epimat) = epinames[rownames(gwmat)]
    emat = epimat # Save non-thresholded var (for diag)
    epimat[epimat > zmax] <- zmax
    epimat[epimat < zmin] <- 0
    epiorder = rownames(epimat)

    # Order
    method = 'ejaccard'
    dt = dist(epimat[,studyord] > zmin, method)
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
    # image(gmat[, studyord], axes=FALSE,col=colred, zlim=c(zmin, zmax))
    # image(epimat[rnmat, studyord], axes=FALSE,col=colred, zlim=c(zmin, zmax))

    # Labels in figure:
    labelsmat = meta[rnmat, 'GROUP']
    faclabels = as.matrix(as.numeric(labelsmat))
    colset = as.character(rdcol$COLOR)
    fix = rnmat[is.na(labelsmat)]
    print(length(fix))
    lablist.mat = label.runs(faclabels, labelsmat, rdcol)

    # Set parameters
    rnorder = rnmat 
    acutnamed = acutmat.nam
    set = tagline 
    lablist = lablist.mat

    png(paste0(imgpref,'top',SHOWGWAS,'_epi.png'),res=450,units='in',width=15,height=17)
    hbreaks = calc.breaks.acut(acutnamed)
    layout(matrix(c(1,2),1,2), widths=c(1.75,10), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    par(mar=c(15, 6, 2, 0.25))
    meta.image(metamat[rnorder,5:1], colvals=colvals, cex=0, horiz=T)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    metaclass = sapply(rev(colnames(metamat)), capitalize)
    text(x=seq(0,1, length.out=ncol(metamat)),
         y=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), 
         labels=metaclass, srt=90, adj=1, xpd=TRUE, cex=.7)
    text(y=lablist[[1]],
         x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
         labels=lablist[[2]], srt=0, adj=1, xpd=TRUE,cex=.8, col=lablist[[3]])
    # Plot gwas alone 
    par(mar=c(15, 0.25, 2, 0.25))
    image(t(epimat[rnorder, studyord]), axes=FALSE,col=colred, zlim=c(zmin, zmax))
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    rows = (seq(0, ncol(gmat), 10) - .5) / (ncol(gmat) - 1)
    abline(h=hbreaks,lty='dotted',lw=1, col='darkgrey')
    abline(v=rows,lty='dotted',lw=1, col='darkgrey')
    text(x=seq(0,1, length.out=length(studyord)),
         y=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), 
         labels=studyord, srt=90, adj=1, xpd=TRUE,cex=.3)
    mtext(tagline, side=3, line=0.25)
    dev.off()


    # Re-order by group
    ordgroup = order(labelsmat, decreasing=TRUE)
    rngroup = rnmat[ordgroup]

    # Params:
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

    # Set parameters
    rnorder = rngroup
    acutnamed = acutgroup.nam
    set = tagline 
    lablist = lablist.group

    png(paste0(imgpref,'top',SHOWGWAS,'_epi_groupord.png'),res=450,units='in',width=15,height=17)
    hbreaks = calc.breaks.acut(acutnamed)
    layout(matrix(c(1,2),1,2), widths=c(1.75,10), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    par(mar=c(15, 6, 2, 0.25))
    meta.image(metamat[rnorder,5:1], colvals=colvals, cex=0, horiz=T)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    metaclass = sapply(rev(colnames(metamat)), capitalize)
    text(x=seq(0,1, length.out=ncol(metamat)),
         y=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), 
         labels=metaclass, srt=90, adj=1, xpd=TRUE, cex=.7)
    text(y=lablist[[1]],
         x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
         labels=lablist[[2]], srt=0, adj=1, xpd=TRUE,cex=.8, col=lablist[[3]])
    # Plot gwas alone 
    par(mar=c(15, 0.25, 2, 0.25))
    image(t(epimat[rnorder, studyord]), axes=FALSE,col=colred, zlim=c(zmin, zmax))
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    rows = (seq(0, ncol(gmat), 10) - .5) / (ncol(gmat) - 1)
    abline(h=hbreaks,lty='dotted',lw=1, col='darkgrey')
    abline(v=rows,lty='dotted',lw=1, col='darkgrey')
    text(x=seq(0,1, length.out=length(studyord)),
         y=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), 
         labels=studyord, srt=90, adj=1, xpd=TRUE,cex=.3)
    mtext(tagline, side=3, line=0.25)
    dev.off()

    # Diagonalize, ignore all <2
    tmat = emat
    tmat[tmat < zmin] = 0
    ll = diag.mat(tmat[rnorder,keep.studies])
    tmp = ll[[1]]
    cto = ll[[3]] # For breaks
    diagord = colnames(tmp)
    vcut =  c(cto[cto ==0], acutnamed[cto])
    vbreaks = calc.breaks.acut(vcut)

    png(paste0(imgpref,'top',SHOWGWAS,'_epi_groupord_diag.png'),res=450,units='in',width=15,height=17)
    hbreaks = calc.breaks.acut(acutnamed)
    layout(matrix(c(1,2),1,2), widths=c(1.75,10), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    par(xpd=FALSE)
    par(mar=c(15, 6, 2, 0.25))
    meta.image(metamat[rnorder,5:1], colvals=colvals, cex=0, horiz=T)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    metaclass = sapply(rev(colnames(metamat)), capitalize)
    text(x=seq(0,1, length.out=ncol(metamat)),
         y=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), 
         labels=metaclass, srt=90, adj=1, xpd=TRUE, cex=.7)
    text(y=lablist[[1]],
         x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
         labels=lablist[[2]], srt=0, adj=1, xpd=TRUE,cex=.8, col=lablist[[3]])
    # Plot gwas alone 
    par(mar=c(15, 0.25, 2, 0.25))
    image(t(epimat[rnorder, diagord]), axes=FALSE,col=colred, zlim=c(zmin, zmax))
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    abline(h=hbreaks,lty='dotted',lw=1, col='darkgrey')
    abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
    text(x=seq(0,1, length.out=length(studyord)),
         y=par()$usr[3]-0.005*(par()$usr[4]-par()$usr[3]), 
         labels=diagord, srt=90, adj=1, xpd=TRUE,cex=.3)
    mtext(tagline, side=3, line=0.25)
    # Add rectangles
    vcind = unique(vcut)[-c(1)] 
    hb = rev(c(par()$usr[3], hbreaks, par()$usr[4]))
    hb = c(hb[vcind], par()$usr[3])
    vb = c(vbreaks, par()$usr[2])
    NV = length(vb)
    rd = cbind(vb[-NV], vb[-1], hb[-NV], hb[-1])
    # One for the full - side:
    rd = rbind(c(par()$usr[1], vb[1], par()$usr[3:4]) ,rd)
    vccols = c('grey',colset[vcind])
    rect(xleft=rd[,1], xright=rd[,2], ybottom=rd[,3], ytop=rd[,4], border=vccols, lwd=1.5)
    par(xpd=NA)
    rect(xleft=rd[,1], xright=rd[,2], ybottom=par()$usr[3], 
         ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
         # border=vccols, col=vccols, lwd=1)
         border='white', col=vccols, lwd=.25)
    dev.off()

    # ------------------------------------
    # Flip and add rectangles on diagonal:
    # take one with strongest overall OR highest PMID number (prob. most recent?) - diff type of enr...
    png(paste0(imgpref,'top',SHOWGWAS,'_epi_groupord_diag_flip.png'),res=450,units='in',width=16,height=17)
    layout(matrix(c(1,2),2,1), heights=c(1.5,10), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    par(xpd=FALSE)
    par(mar=c(0, 15, 6, 0.25))
    meta.image(metamat[rnorder,], colvals=colvals, cex=0, horiz=F)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    metaclass = sapply(colnames(metamat), capitalize)
    text(y=seq(0,1, length.out=ncol(metamat)),
         x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
         labels=metaclass, srt=0, adj=1, xpd=TRUE, cex=.7)
    text(x=lablist[[1]],
         y=par()$usr[4]+0.05*(par()$usr[4]-par()$usr[3]), 
         labels=lablist[[2]], srt=90, adj=0, xpd=TRUE,cex=.8, col=lablist[[3]])
    # Plot gwas alone 
    par(mar=c(.25, 15, .25, 0.25))
    image(epimat[rnorder, rev(diagord)], axes=FALSE,col=colred, zlim=c(zmin, zmax))
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    hbreaks = calc.breaks.acut(acutnamed)
    vbreaks = calc.breaks.acut(vcut)
    abline(v=hbreaks,lty='dotted',lw=1, col='darkgrey')
    abline(h=par()$usr[2] - vbreaks,lty='dotted',lw=1, col='darkgrey')
    text(y=seq(0,1, length.out=length(diagord)),
         x=par()$usr[1]-0.005*(par()$usr[2]-par()$usr[1]), 
         labels=rev(diagord), srt=0, adj=1, xpd=TRUE,cex=.3)
    # mtext(tagline, side=3, line=0.25)
    # Add rectangles
    vcind = unique(vcut)[-c(1)] 
    hb = rev(c(par()$usr[1], hbreaks, par()$usr[2]))
    hb = c(hb[vcind], par()$usr[1])
    vb = c(vbreaks, par()$usr[4])
    NV = length(vb)
    rd = cbind(vb[-NV], vb[-1], hb[-NV], hb[-1])
    # One for the full - side:
    rd = rbind(c(par()$usr[1], vb[1], par()$usr[3:4]) ,rd)
    vccols = c('grey',colset[vcind])
    rect(xleft=rd[,3], xright=rd[,4], ybottom=par()$usr[2] - rd[,1], 
         ytop=par()$usr[2] - rd[,2], border=vccols, lwd=1)
    par(xpd=NA)
    rect(ybottom=par()$usr[2] - rd[,1], ytop=par()$usr[2] - rd[,2],
         xleft=par()$usr[1], xright=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
         border='white', col=vccols, lwd=.25)
    dev.off()

}

