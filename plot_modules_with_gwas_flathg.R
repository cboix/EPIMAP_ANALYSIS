#!/usr/bin/R
# -----------------------------------------------------
# NOTE: Same as plot_modules_with_gwas but for the 
# flat module enrichments.
# 
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
source(paste0(bindir, 'great_enrichment/auxfunctions_GREAT.R'))
source(paste0(bindir, 'load_metadata.R'))
today <- format(Sys.time(), "%m%d%y")
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

# ---------------------
# Defaults for modules:
# ---------------------
filepref = 'cls_merge2_wH3K27ac100_300'
tagline = 'ChromHMM Enhancer Modules'
imgdir = paste0(img, "clusters/") 

# Load in and process data:
commandArgs <- function(trailingOnly=TRUE){
    c(filepref, tagline, imgdir) }
source(paste0(bindir, 'load_modules_data.R'))

# -----------------------
# Defaults for gwas data:
# -----------------------
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = FALSE
use.adj = TRUE
use.strict = FALSE
# use.strict = TRUE

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only, use.adj, use.strict) }
source(paste0(bindir, 'load_validation_gwastree_enrichments.R'))

# Prefix:
imgpref = paste0(img, 'clusters/clusters_')

# --------------------------------
# Add gwas enrichment information:
# --------------------------------
# Filter down to only >20k ind.
NIND = 20000 # Remove all the really low GWAS
keptgw = gwssdf[gwssdf$sampsize > NIND,]
rmat = regmats[['mod']]
rmat[is.na(rmat)] = 0
cutoff = 3
kuid = names(which(apply(rmat > cutoff, 1, sum) > 0))
kuid = kuid[kuid %in% keptgw$uid]
gwmat = t(rmat[kuid,])
rownames(gwmat) = paste0('c', c(1:nrow(gwmat) - 1))
gwmat[gwmat < 1] = 0
keep.cls = list()
keep.cls[[tagline]] = rownames(gwmat)

# Choose # gwas to show:
# SHOWGWAS=300
SHOWGWAS=length(kuid)
zmax = 20
zmin=3
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

plot.gwas <-  function(gwas, set, ordll, zmin=1, zmax=5, ablwd=1,
                       gwasord=NULL, labeled=TRUE, calc.only=FALSE){
    # GWAS PLOT HERE:
    ord = ordll[[1]]
    vbreaks = ordll[[2]]
    clscuts = ordll[[3]]
    subord = ord[ord %in% rownames(gwas)]
    gmat = gwas[subord,]
    if (!is.null(gwasord)){ gmat = gmat[,gwasord] }
    # Diagonalize, ignoring lower than cutoff
    gmat[gmat < zmin] <- 0
    ll = diag.mat2(gmat)
    gmat = ll[[1]]
    cto = ll[[3]] # For breaks
    vcut =  c(cto[cto ==0], clscuts[cto])
    hbreaks = calc.breaks.acut(vcut)
    # Threshold for plotting:
    gmat[gmat > zmax] <- zmax
    if (!calc.only){
        image(gmat, axes=FALSE,col=colred, zlim=c(zmin, zmax), useRaster=TRUE)
        abline(h=par()$usr[3:4],lty=1,lw=0.5)
        abline(v=par()$usr[1:2],lty=1,lw=0.5)
        abline(v=vbreaks,lty='dotted',lw=ablwd, col='darkgrey')
        abline(h=hbreaks,lty='dotted',lw=ablwd, col='darkgrey')
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
png(paste0(imgpref, "gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord.png"),res=450,units='in',width=12,height=17)
ratio = 1
layout(matrix(c(1:12),4,3), heights=c(.5,8,8 * ratio, 1), widths=c(1.25,8,.5), TRUE)
par(xpd=FALSE)
par(yaxs="i")
par(xaxs="i")
dist.order = list()
# Metadata matrix:
par(mar=c(0.25, 6, 0, 0))
image(t(as.matrix(metamat) == ""), col='white', axes=F, useRaster=TRUE)
metaclass = sapply(rev(colnames(metamat)), capitalize)
text(x=seq(0,1, length.out=ncol(metamat)),
     y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), 
     labels=metaclass, srt=90, adj=0, xpd=TRUE, cex=.7)
par(mar=c(.25, 6, 0, 0))
meta.image(metamat[rnorder,5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
text(y=lablist[[1]],
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lablist[[2]], srt=0, adj=1, xpd=TRUE,cex=.8, col=lablist[[3]])
abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
# GWAS labels:
par(mar=c(0.25, 6, 0, 0))
par(xpd=NA)
image(gwmat[,gwnam], col='white', axes=F, useRaster=TRUE)
text(y=seq(0,1, length.out=length(gwnam)),
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=gwnam, srt=0, adj=1, xpd=TRUE,cex=.25)
par(xpd=FALSE)
# Add labels
plot.new()
par(mar=c(0.25, 0.25, 2, 0.25))
dist.order = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, calc.only=TRUE, subset=TRUE)
clsord = dist.order[[1]]
meta.image(enrichmat[clsord,5:1], colvals=colvals, cex=0, horiz=F, useRaster=TRUE)
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
xloc = -0.16
rect(ybottom=rectdf$y1+0.0001, ytop=rectdf$y2-0.0001,
     xright=par()$usr[1], xleft=xloc, 
     border=vccols, col=NA, lwd=.25)
rectdf$ym = (rectdf$y1 + rectdf$y2) / 2
rectdf$COLOR = vccols
rectdf = merge(rectdf, rdcol, all.x=TRUE)
rectdf$GROUP = as.character(rectdf$GROUP)
rectdf$GROUP[is.na(rectdf$GROUP)] = 'None'
text(x=xloc, y=rectdf$ym, labels=rectdf$GROUP, 
     col=rectdf$COLOR, cex=.75, adj=1)
par(xpd=FALSE)
par(mar=c(0.25, 0.25, 0, 0.25))
plot.counts(counts, set, dist.order[[set]])
plot.new()
# Availability:
par(mar=c(.25, 0, 0, 0.25))
avail = as.matrix(wm[main.marks, rnorder])
image(avail, axes=F, col=c('white', 'darkgrey'), useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
par(mar=c(.25, 0, 0, 0.25))
image(avail, axes=F, col='white', useRaster=TRUE)
text(x=seq(0,1, length.out=length(main.marks)),
     y=par()$usr[4], 
     labels=main.marks, srt=90, adj=1, xpd=TRUE,cex=.5)
dev.off()



# --------------------------------------
# Plot only the GWAS results on modules:
# --------------------------------------
height = 11 * ncol(finalmat) / 500
termcex = 0.18


# png(paste0(imgpref, "gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_onlygwas.png"),res=450,units='in',width=8,height=height)
pdf(paste0(imgpref, "gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_onlygwas.pdf"),width=8,height=height)
ratio = 1
layout(matrix(c(1:2),1,2), widths=c(4,7), TRUE)
par(xpd=FALSE)
par(yaxs="i")
par(xaxs="i")
dist.order = list()
sp = 0.15
# GWAS labels:
par(mar=c(sp, 6, sp, sp))
par(xpd=NA)
image(gwmat[,gwnam], col='white', axes=F, useRaster=TRUE)
text(y=seq(0,1, length.out=length(gwnam)),
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
     # labels=gwnam, srt=0, adj=1, xpd=TRUE,cex=termcex)
     labels=gwnam, srt=0, adj=1, xpd=TRUE,cex=.2)
par(xpd=FALSE)
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
gll = plot.gwas(finalmat, set, dist.order[[set]], gwasord=studyord, ablwd=.5)
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


# ---------------------
# Split into two plots:
# ---------------------
plot_allgwas_splitid = function(pltid, sp=.1, txtcex=.15){
    ratio = 1
    layout(matrix(c(1:2),1,2), widths=c(4,7), TRUE)
    par(xpd=FALSE)
    par(yaxs="i")
    par(xaxs="i")
    dist.order = list()
    # GWAS labels:
    par(mar=c(sp, 6, sp, sp))
    par(xpd=NA)
    image(gwmat[,gwnam[pltid]], col='white', axes=F, useRaster=TRUE)
    gcol = c('grey', colset)[gcut + 1]
    text(y=seq(0,1, length.out=length(gwnam[pltid])),
         x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
         labels=gwnam[pltid], srt=0, adj=1, xpd=TRUE,
         cex=txtcex, col=gcol[pltid])
    par(xpd=FALSE)
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
    # gll = plot.gwas(finalmat, set, dist.order[[set]], gwasord=studyord, ablwd=.5)
    # GWAS PLOT (for split plot):
    ordll = dist.order[[set]]
    ord = ordll[[1]]
    vbreaks = ordll[[2]]
    clscuts = ordll[[3]]
    gwas = gwmat
    subord = ord[ord %in% rownames(gwas)]
    gmat = gwas[subord,studyord]
    # Diagonalize, ignoring lower than cutoff
    gmat[gmat < zmin] <- 0
    ll = diag.mat2(gmat)
    gmat = ll[[1]]
    cto = ll[[3]][pltid] # For breaks
    vcut =  c(cto[cto ==0], clscuts[cto])
    hbreaks = calc.breaks.acut(vcut)
    # Threshold for plotting:
    zmin=1
    zmax=5
    gmat[gmat > zmax] <- zmax
    image(gmat[,pltid], axes=FALSE,col=colred, zlim=c(zmin, zmax), useRaster=TRUE)
    ablwd=0.5
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    abline(v=vbreaks,lty='dotted',lw=ablwd, col='darkgrey')
    abline(h=hbreaks,lty='dotted',lw=ablwd, col='darkgrey')
    rll = calc.breaks.rect(hcls=gcut[pltid], vcls=dist.order[[set]][[3]], colset)
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
}

midpt = round(ncol(gwmat)/2)
botid = 1:midpt
topid = (midpt + 1):ncol(gwmat)

pdf(paste0(imgpref, "gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_onlygwas_tophalf.pdf"), width=8.25, height=10)
plot_allgwas_splitid(topid)
dev.off()

pdf(paste0(imgpref, "gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_onlygwas_bothalf.pdf"), width=8.25, height=10)
plot_allgwas_splitid(botid)
dev.off()



# --------------------------------------------------------------
# Reduce and show only a couple of GWAS (representative - wdlt):
# --------------------------------------------------------------
# Include best gwas signal.
NTERMS = 200 # total number terms
termcex = 0.35
# NTERMS = 50 # total number terms
# termcex = 0.4
vcut = gcut
NREP = NTERMS - length(unique(vcut)) # 1 per at least.
rdf = aggregate(count ~ cls, data.frame(cls = factor(vcut, levels =unique(vcut)), count=1), sum)
# Add up to NTERMS representatively:
rrp = rdf$count / (length(vcut) / NREP)
rdf$rep = floor(rrp) + 1
add1 = head(order(rrp %% 1, decreasing=T), NTERMS - sum(rdf$rep))
rdf$rep[add1] = rdf$rep[add1] + 1
rdf$cls = as.numeric(as.character(rdf$cls))
diagord = gwnam
diagord2 = diagord
diagord2[nchar(diagord2) > 150] = ''
term_words <- strsplit(diagord2, "[ _,.]");
tab_all <- sort(table(unlist(term_words)));
terms = c()
wcls = c()
for (i in 1:nrow(rdf)){
    wind = which(vcut == rdf$cls[i])
    subrn = diagord2[wind]
    # if (length(subrn) > 5){ filt=TRUE } else { filt=FALSE }
    if (length(subrn) > 3){ filt=TRUE } else { filt=FALSE }
    bs = round(length(subrn) / rdf$rep[i])
    if (bs == 0){bs = 1}
    sterms = get_summary_terms(subrn, binsize=bs, filtering=filt, tab_all=tab_all)
    if (rdf$rep[i] != length(sterms)){
        print(i)
    }
    terms = c(terms, sterms)
    wcls = c(wcls, rep(rdf$cls[i], length(sterms)))
}

# Color wordlets by category:
wcol = c('grey', colset)[wcls + 1]
wordlets = terms
wind = sapply(wordlets, function(x){which(diagord2 == x)})

# ------------------------------------
# Flip and add rectangles on diagonal:
# png(paste0(imgpref, "gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_repr.png"),res=450,units='in',width=7.75,height=5)
pdf(paste0(imgpref, "gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_repr.pdf"),width=7.75,height=10)
sp = 0.15
layout(matrix(c(1,2),1,2), widths=c(4.5,7), TRUE)
par(yaxs="i")
par(xaxs="i")
par(xpd=FALSE)
# Plot labels:
par(mar=c(sp, 6, sp, sp))
image(gwmat[,gwnam], col='white', axes=F, useRaster=TRUE)
par(xpd=TRUE)
rx=c(0.002,0.045, 0.2,0.25, 0.255) * .3
xl = wordlets
xc = wcol
xx = seq(0+1e-3, 1-1e-3, length.out=length(xl))
xw = wind / (length(diagord) +1)
x = par()$usr[2]-rx*(diff(par()$usr[1:2]))
text(y=xx, x=x[5], labels=xl, srt=0, adj=1, xpd=TRUE, cex=termcex, col=wcol)
segments(x0=x[1], y0=xw, x1=x[2], y1=xw, col=wcol, lwd=.5)
segments(x0=x[2], y0=xw, x1=x[3], y1=xx, col=wcol, lwd=.5)
segments(x0=x[3], y0=xx, x1=x[4], y1=xx, col=wcol, lwd=.5)
par(xpd=FALSE)
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
gll = plot.gwas(finalmat, set, dist.order[[set]], gwasord=studyord, ablwd=.5)
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

