#!/usr/bin/R
# -----------------------------------------------------
# Plot modules from CHMM/other clustering with GWAS
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
gwasfile = 'observed_aux_18_on_mixed_impobs_QCUT_PROM_bin_on_mixed_impobs_5000_enrich.tsv'
# gwasfile = 'observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_5000_enrich.tsv'
# gwasfile = 'observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_1000_enrich.tsv'
# gwasfile = 'observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_0_enrich.tsv'
# filepref = 'cls_merge2_wH3K27ac100_300'
#tagline = 'ChromHMM Enhancers'
filepref = 'prom_cls/cls_merge2_wH3K27ac100_300'
tagline = 'ChromHMM Promoters'
extension = 5000
imgdir = paste0(img, "clusters/") 

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need clusters filename")
} else {        
    gwasfile = args[1]
    filepref = args[2]
    tagline = args[3]
    extension = args[4]
    if (length(args) > 4){ 
        imgdir = args[5] 
    }
}

# -------------------------
# Load in and process data:
# -------------------------
commandArgs <- function(trailingOnly=TRUE){
    c(filepref, tagline, imgdir) }
source(paste0(bindir, 'load_modules_data.R'))

# Prefix:
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
# imgpref = paste0(imgdir, 'clusters_gwas_', filepref, '_e', extension, '_')
imgpref = paste0(imgdir, 'clusters_gwas_', sub("/","_", filepref), '_e', extension, '_')

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
keep.cls = list()
keep.cls[[tagline]] = rownames(gwmat)

# Choose # gwas to show:
SHOWGWAS=300
zmax = 5
zmin=0.5
# gwasmarg = sort(apply(gwmat, 2, sum), decreasing=T)
gwasmarg = sort(apply(gwmat, 2, max), decreasing=T)
keep.studies = names(head(gwasmarg, SHOWGWAS))
# Order the top studies:
r2 = reord(t(gwmat[, keep.studies]) > zmin, 'Jaccard')
studyord = rownames(r2)
r3 = reord(gwmat[, keep.studies] > zmin, 'eJaccard')
roword = rownames(r3)
finalmat = gwmat[roword, studyord]

# Threshold for plotting:
# Plot diagonalized -> first order by cluster (rows) -> diag
gmat = finalmat

plot.gwas <-  function(gwas, set, ordll, zmin=1, zmax=5, gwasord=NULL, labeled=TRUE, calc.only=FALSE){
    # GWAS PLOT HERE:
    ord = ordll[[1]]
    vbreaks = ordll[[2]]
    clscuts = ordll[[3]]
    subord = ord[ord %in% rownames(gwas)]
    gmat = gwas[subord,]
    if (!is.null(gwasord)){ gmat = gmat[,gwasord] }
    # Diagonalize, ignoring lower than cutoff
    gmat[gmat < zmin] <- 0
    ll = diag.mat(gmat)
    gmat = ll[[1]]
    cto = ll[[3]] # For breaks
    vcut =  c(cto[cto ==0], clscuts[cto])
    hbreaks = calc.breaks.acut(vcut)
    # Threshold for plotting:
    gmat[gmat > zmax] <- zmax
    if (!calc.only){
        image(gmat, axes=FALSE,col=colred, zlim=c(zmin, zmax))
        abline(h=par()$usr[3:4],lty=1,lw=0.5)
        abline(v=par()$usr[1:2],lty=1,lw=0.5)
        abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
        abline(h=hbreaks,lty='dotted',lw=1, col='darkgrey')
    }
    return(list(colnames(gmat), vcut))
}

# Set groupdist parameters
rnorder = rngroup # ord.gwas
acutnamed = acutgroup.nam
set = tagline 
lablist = lablist.group

dist.order = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, subset=TRUE, calc.only=TRUE)
gll = plot.gwas(finalmat, set, dist.order, gwasord=studyord, calc.only=TRUE)
gwnam = gll[[1]]
gcut = gll[[2]]


# Make plot:
png(paste0(imgpref,'top', SHOWGWAS,'_groupord.png'),res=450,units='in',width=12,height=17)
ratio = 1
layout(matrix(c(1:12),4,3), heights=c(.5,8,8 * ratio, 1), widths=c(1.25,8,.5), TRUE)
par(xpd=FALSE)
par(yaxs="i")
par(xaxs="i")
dist.order = list()
# Metadata matrix:
par(mar=c(0.25, 6, 0, 0))
image(t(as.matrix(metamat) == ""), col='white', axes=F)
metaclass = sapply(rev(colnames(metamat)), capitalize)
text(x=seq(0,1, length.out=ncol(metamat)),
     y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), 
     labels=metaclass, srt=90, adj=0, xpd=TRUE, cex=.7)
par(mar=c(.25, 6, 0, 0))
meta.image(metamat[rnorder,5:1], colvals=colvals, cex=0, horiz=T)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
text(y=lablist[[1]],
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lablist[[2]], srt=0, adj=1, xpd=TRUE,cex=.8, col=lablist[[3]])
abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
# GWAS labels:
par(mar=c(0.25, 6, 0, 0))
par(xpd=NA)
image(gwmat[,gwnam], col='white', axes=F)
text(y=seq(0,1, length.out=length(gwnam)),
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=gwnam, srt=0, adj=1, xpd=TRUE,cex=.25)
par(xpd=FALSE)
# Add labels
plot.new()
par(mar=c(0.25, 0.25, 2, 0.25))
dist.order = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, calc.only=TRUE, subset=TRUE)
clsord = dist.order[[1]]
meta.image(enrichmat[clsord,5:1], colvals=colvals, cex=0, horiz=F)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
vbreaks = dist.order[[set]][[2]]
abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
mtext(set, side=3, cex=1.3)
# Plot clusters and counts
set = tagline
par(mar=c(0.25, 0.25, 0, 0.25))
dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, title=FALSE, subset=TRUE)
abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
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
# GWAS plot
par(mar=c(0.25, 0.25, 0, 0.25))
gll = plot.gwas(finalmat, set, dist.order[[set]], gwasord=studyord)
gwnam = gll[[1]]
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
     border=vccols, lwd=1)
par(xpd=NA)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
rect(xleft=rectdf$x1, xright=rectdf$x2, ybottom=par()$usr[3], 
     ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=vccols, lwd=.25)
par(xpd=FALSE)
par(mar=c(0.25, 0.25, 0, 0.25))
plot.counts(counts, set, dist.order[[set]])
plot.new()
# Availability:
par(mar=c(.25, 0, 0, 0.25))
avail = as.matrix(wm[main.marks, rnorder])
image(avail, axes=F, col=c('white', 'darkgrey'))
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
par(mar=c(.25, 0, 0, 0.25))
image(avail, axes=F, col='white')
text(x=seq(0,1, length.out=length(main.marks)),
     y=par()$usr[4], 
     labels=main.marks, srt=90, adj=1, xpd=TRUE,cex=.5)
dev.off()


