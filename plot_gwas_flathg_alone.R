#!/usr/bin/R
# -----------------------------------------------------
# Plot epigenomes (or modules) with GWAS
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
library(UpSetR)

# Defaults:
# gwasfile = 'observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_5000_enrich.tsv'
# extension = 5000
filepref = 'cls_merge2_wH3K27ac100_raw'
tagline = 'ChromHMM Enhancers Epigenomes'
imgdir = paste0(img, "clusters/") 

# -----------------------
# Defaults for gwas data:
# -----------------------
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = FALSE
use.adj = TRUE
use.strict = FALSE

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only, use.adj, use.strict) }
source(paste0(bindir, 'load_validation_gwastree_enrichments.R'))

# Plot prefix:
imgpref = paste0(img, 'clusters/')

# --------------------------------
# Add gwas enrichment information:
# --------------------------------
# Filter down to only >20k ind.
NIND = 20000 # Remove all the really low GWAS
keptgw = gwssdf[gwssdf$sampsize >= NIND,]
rmat = regmats[['epi']]  # Epigenomes flat matrix
rmat[is.na(rmat)] = 0
kuid = names(which(apply(rmat > 0, 1, sum) > 0))
kuid = kuid[kuid %in% keptgw$uid]
gwmat = t(rmat[kuid,])
rownames(gwmat) = paste0('c', c(1:nrow(gwmat) - 1))
gwmat[gwmat < 1] = 0

dim(gwdf[gwdf$uid %in% keptgw$uid,])

# Choose # gwas to show:
# SHOWGWAS=350
SHOWGWAS=800 # Cover all.
# gwasmarg = sort(apply(gwmat, 2, sum), decreasing=T)
gwasmarg = sort(apply(gwmat, 2, max), decreasing=T)
keep.studies = names(head(gwasmarg, SHOWGWAS))
zmax = 20
zmin=3
print(sum(gwmat[,keep.studies] > zmin))

# Order the top studies:
r2 = reord(t(gwmat[, keep.studies]) > zmin, 'Jaccard')
studyord = rownames(r2)
r3 = reord(gwmat[, keep.studies] > zmin, 'eJaccard')
roword = rownames(r3)
finalmat = gwmat[roword, studyord]
namesfile = paste0(filepref, '_names.tsv')
epinames = scan(namesfile,'c')
clsn = paste0('c', 1:length(epinames) - 1)
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


# Establish levels/ordering:
odf = rdcol  # Copy don't change rdcol
odf = rbind(odf, c('None','grey80','Other'))
clvs = c('Tissue', 'Primary Cell', 'Immune', 'Neural', 'Stem-like', 'Other')
odf$category = factor(odf$category, levels=clvs)
# Sort alpha then group:
odf = odf[order(odf$GROUP),]
odf = odf[order(odf$category),]

# Labels in figure:
labelsmat = meta[rnmat, 'GROUP']
labelsmat = as.character(labelsmat)
labelsmat = factor(labelsmat, levels=odf$GROUP)
faclabels = as.matrix(as.numeric(labelsmat))
colset = as.character(rdcol$COLOR)
fix = rnmat[is.na(labelsmat)]
print(length(fix))
lablist.mat = label.runs(faclabels, labelsmat, rdcol)

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

# Diagonalize, ignore all <2
tmat = emat
tmat[tmat < zmin] = 0
ll = diag.mat2(tmat[rnorder,keep.studies], ratio=.25)
tmp = ll[[1]]
cto = ll[[3]] # For breaks
diagord = colnames(tmp)
vcut =  c(cto[cto ==0], acutnamed[cto])
vbreaks = calc.breaks.acut(vcut)
full.diagord = diagord
full.vcut = vcut

png(paste0(imgpref, "epigenomes_gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord.png"),res=450,units='in',width=17,height=17)
hbreaks = calc.breaks.acut(acutnamed)
layout(matrix(c(1,2),1,2), widths=c(1.75,10), TRUE)
par(yaxs="i")
par(xaxs="i")
par(xpd=FALSE)
par(mar=c(15, 6, 2, 0.25))
meta.image(metamat[rnorder,5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
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
image(t(epimat[rnorder, diagord]), axes=FALSE,col=colred, zlim=c(zmin, zmax), useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
abline(h=hbreaks,lty='dotted',lw=1, col='darkgrey')
abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
text(x=seq(0,1, length.out=length(studyord)),
     y=par()$usr[3]-0.005*(par()$usr[4]-par()$usr[3]), 
     labels=diagord, srt=90, adj=1, xpd=TRUE,cex=.3)
mtext(tagline, side=3, line=0.25)
# Add rectangles
rll = calc.breaks.rect(hcls=acutnamed, vcls=vcut, colset)
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
dev.off()


# ----------------------------
# Flip figure and diagonalize:
# ----------------------------
# png(paste0(imgpref, "epigenomes_gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_diag_flip.png"),res=450,units='in',width=12,height=17)
pdf(paste0(imgpref, "epigenomes_gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_diag_flip.pdf"), width=12, height=17)
layout(matrix(c(1,2),2,1), heights=c(1.2,10), TRUE)
par(yaxs="i")
par(xaxs="i")
par(xpd=FALSE)
par(mar=c(0, 15, 7, 0.25))
meta.image(metamat[rnorder,], colvals=colvals, cex=0, horiz=F, useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
metaclass = sapply(colnames(metamat), capitalize)
text(y=seq(0,1, length.out=ncol(metamat)),
     x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=metaclass, srt=0, adj=1, xpd=TRUE, cex=.6)
# Plot labels:
box.pad = 0.02
xw = lablist[[1]]
xx = space.1d(xw, box.pad=box.pad)
xx = space.1d(xx, box.pad=box.pad)
rx = c(0.05, 0.15, 0.4, 0.45)
x = par()$usr[4]+rx*(diff(par()$usr[3:4]))
text(x=xx, y=x[4], labels=lablist[[2]], col=lablist[[3]],
     srt=90, adj=0, xpd=TRUE, cex=.8)
par(xpd=TRUE)
segments(x0=xw, y0=x[2], x1=xx, y1=x[3], col=lablist[[3]])
segments(x0=xw, y0=x[1], x1=xw, y1=x[2], col=lablist[[3]])
par(xpd=FALSE)
par(mar=c(.25, 15, .25, 0.25))
image(epimat[rnorder, rev(diagord)], axes=FALSE,col=colred, zlim=c(zmin, zmax), useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
text(y=seq(0,1, length.out=length(diagord)),
     x=par()$usr[1]-0.005*(par()$usr[2]-par()$usr[1]), 
     labels=rev(diagord), srt=0, adj=1, xpd=TRUE,cex=.2)
# Add rectangles
rll = calc.breaks.rect(hcls=rev(vcut), vcls=acutnamed, colset)
rectdf = rll[[1]]
vccols = rll[[2]]
rord = order(rectdf$y1)
rectdf = rectdf[rord,]
vccols = vccols[rord]
hbreaks = calc.breaks.acut(acutnamed)
abline(h=unique(c(rectdf$y1, rectdf$y2)),lty='dotted', lwd=.5, col='darkgrey')
abline(v=hbreaks,lty='dotted', lwd=.5, col='darkgrey')
# Add grey - ubq rectangle:
rectdf = rbind(c(x1=par()$usr[1], x2=par()$usr[2],
                 y1=rectdf$y2[nrow(rectdf)], y2=par()$usr[4]), rectdf)
vccols = c('grey',vccols)
rect(xleft=rectdf$x1, xright=rectdf$x2,
     ybottom=rectdf$y1, ytop=rectdf$y2,
     border=vccols, lwd=1)
par(xpd=NA)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
dev.off()


# --------------------------------------------------------------
# Reduce and show only a couple of GWAS (representative - wdlt):
# --------------------------------------------------------------
# Include best gwas signal.
NTERMS = 100 # total number terms
NTERMS = 107 # total number terms
termcex = 0.3
# NTERMS = 50 # total number terms
# termcex = 0.4
NREP = NTERMS - length(unique(vcut)) # 1 per at least.
rdf = aggregate(count ~ cls, data.frame(cls = factor(vcut, levels =unique(vcut)), count=1), sum)
# Add up to NTERMS representatively:
rrp = rdf$count / (length(vcut) / NREP)
rdf$rep = floor(rrp) + 1
add1 = head(order(rrp %% 1, decreasing=T), NTERMS - sum(rdf$rep))
rdf$rep[add1] = rdf$rep[add1] + 1
rdf$cls = as.numeric(as.character(rdf$cls))
diagord2 = diagord
diagord2[sapply(diagord2, nchar) > 100] = ''
term_words <- strsplit(diagord2, "[ _,.]");
tab_all <- sort(table(unlist(term_words)));
terms = c()
wcls = c()
for (i in 1:nrow(rdf)){
    wind = which(vcut == rdf$cls[i])
    subrn = diagord2[wind]
    # if (length(subrn) > 5){ filt=TRUE } else { filt=FALSE }
    if (length(subrn) > 3){ filt=TRUE } else { filt=FALSE }
    sterms = get_summary_terms(subrn, binsize=round(length(subrn) / rdf$rep[i]), 
                               filtering=filt, tab_all=tab_all)
    terms = c(terms, sterms)
    wcls = c(wcls, rep(rdf$cls[i], length(sterms)))
}
# Color wordlets by category:
wcol = c('grey', colset)[rev(wcls + 1)]
wordlets = rev(terms)
wind = unlist(sapply(wordlets, function(x){which(rev(diagord2) == x)}))
length(wind)


gw.col_fun = function(x, palette=colred){
    bin <- cut(x, seq(zmin, zmax, length.out=length(palette)), include.lowest=T) 
    palette[bin] }
gw.legend = Legend(at = c(3, 10, 15, 20), 
                   labels_gp = gpar(fontsize=5),
                   title_gp = gpar(fontsize=5, fontface='bold'),
                   col_fun=gw.col_fun, title_position = "topleft", 
                   title="-log10p", direction = 'vertical')
# colred, zlim=c(zmin, zmax)
plegend = packLegend(gw.legend)


# png(paste0(imgpref, "epigenomes_gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_diag_flip_repr.png"),res=450,units='in',width=7.75,height=5)
pdf(paste0(imgpref, "epigenomes_gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_diag_flip_repr.pdf"), width=8, height=5)
sp = 0.15
rsp = 1.25
lsp = 13.5
layout(matrix(c(1,2),2,1), heights=c(2,8), TRUE)
par(yaxs="i")
par(xaxs="i")
par(xpd=FALSE)
par(mar=c(0, lsp, 3.75, rsp))
meta.image(metamat[rnorder,], colvals=colvals, cex=0, horiz=F, useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
metaclass = sapply(colnames(metamat), capitalize)
text(y=seq(0,1, length.out=ncol(metamat)),
     x=par()$usr[2]+0.001*(par()$usr[2]-par()$usr[1]), 
     labels=metaclass, srt=0, adj=0, xpd=TRUE, cex=.3)
# Plot labels:
box.pad = 0.02
xw = lablist[[1]]
xx = space.1d(xw, box.pad=box.pad)
xx = space.1d(xx, box.pad=box.pad)
rx = c(0.05, 0.15, 0.4, 0.45)
x = par()$usr[4]+rx*(diff(par()$usr[3:4]))
text(x=xx, y=x[4], labels=lablist[[2]], col=lablist[[3]],
     srt=90, adj=0, xpd=TRUE, cex=.4)
par(xpd=TRUE)
segments(x0=xw, y0=x[2], x1=xx, y1=x[3], col=lablist[[3]])
segments(x0=xw, y0=x[1], x1=xw, y1=x[2], col=lablist[[3]])
par(xpd=FALSE)
par(mar=c(sp, lsp, sp, rsp))
image(epimat[rnorder, rev(diagord)], axes=FALSE,col=colred, zlim=c(zmin, zmax), useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
# add wordlet text:
# par(xpd=TRUE)
par(xpd=NA)
rx=c(0.005,0.01, 0.03, 0.034, 0.035)
x = par()$usr[1]-rx*(diff(par()$usr[1:2]))
xl = wordlets
xc = wcol
xx = seq(0 + 1e-3, 1.25 - 1e-3, length.out=length(xl))
xw = wind / (length(diagord) +1)
text(y=xx, x=x[5], labels=xl,
     srt=0, adj=1, xpd=NA, cex=termcex, col=wcol)
segments(x0=x[1], y0=xw, x1=x[2], y1=xw, col=wcol, lwd=.5)
segments(x0=x[2], y0=xw, x1=x[3], y1=xx, col=wcol, lwd=.5)
segments(x0=x[3], y0=xx, x1=x[4], y1=xx, col=wcol, lwd=.5)
par(xpd=FALSE)
# Add rectangles
rll = calc.breaks.rect(hcls=rev(vcut), vcls=acutnamed, colset)
rectdf = rll[[1]]
vccols = rll[[2]]
rord = order(rectdf$y1)
rectdf = rectdf[rord,]
vccols = vccols[rord]
hbreaks = calc.breaks.acut(acutnamed)
abline(h=unique(c(rectdf$y1, rectdf$y2)),lty='dotted', lwd=.5, col='darkgrey')
abline(v=hbreaks,lty='dotted', lwd=.5, col='darkgrey')
# Add grey - ubq rectangle:
rectdf = rbind(c(x1=par()$usr[1], x2=par()$usr[2],
                 y1=rectdf$y2[nrow(rectdf)], y2=par()$usr[4]), rectdf)
vccols = c('grey',vccols)
rect(xleft=rectdf$x1, xright=rectdf$x2,
     ybottom=rectdf$y1, ytop=rectdf$y2,
     border=vccols, lwd=.75)
par(xpd=NA)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
draw(plegend, x = unit(0.25,'in'), y=unit(4.75,'in'), just = "top")
dev.off()


# --------------------------------------------
# Reduce the GWAS matrix to top/large studies:
# --------------------------------------------
NIND = 50000 # Remove all the really low GWAS
keptgw = gwssdf[gwssdf$sampsize > NIND,]
kuid = colnames(emat)
kuid = kuid[kuid %in% keptgw$uid]

# Diagonalize, ignore all lt 3
tmat = emat
tmat[tmat < zmin] = 0
ll = diag.mat2(tmat[rnorder, kuid], ratio=.25)
tmp = ll[[1]]
cto = ll[[3]] # For breaks
diagord = colnames(tmp)
vcut =  c(cto[cto ==0], acutnamed[cto])
vbreaks = calc.breaks.acut(vcut)

# Include best gwas signal.
NTERMS = 100 # total number terms
termcex = 0.3
# NTERMS = 50 # total number terms
# termcex = 0.4
NREP = NTERMS - length(unique(vcut)) # 1 per at least.
rdf = aggregate(count ~ cls, data.frame(cls = factor(vcut, levels =unique(vcut)), count=1), sum)
# Add up to NTERMS representatively:
rrp = rdf$count / (length(vcut) / NREP)
rdf$rep = floor(rrp) + 1
add1 = head(order(rrp %% 1, decreasing=T), NTERMS - sum(rdf$rep))
rdf$rep[add1] = rdf$rep[add1] + 1
rdf$cls = as.numeric(as.character(rdf$cls))
diagord2 = diagord
diagord2[sapply(diagord2, nchar) > 100] = ''
term_words <- strsplit(diagord2, "[ _,.]");
tab_all <- sort(table(unlist(term_words)));
terms = c()
wcls = c()
for (i in 1:nrow(rdf)){
    wind = which(vcut == rdf$cls[i])
    subrn = diagord2[wind]
    # if (length(subrn) > 5){ filt=TRUE } else { filt=FALSE }
    if (length(subrn) > 3){ filt=TRUE } else { filt=FALSE }
    sterms = get_summary_terms(subrn, binsize=round(length(subrn) / rdf$rep[i]), 
                               filtering=filt, tab_all=tab_all)
    terms = c(terms, sterms)
    wcls = c(wcls, rep(rdf$cls[i], length(sterms)))
}

# Color wordlets by category:
wcol = c('grey', colset)[rev(wcls + 1)]
wordlets = rev(terms)
wind = sapply(wordlets, function(x){which(rev(diagord2) == x)})

# ------------------------------------
# Flip and add rectangles on diagonal:
png(paste0(imgpref, "epigenomes_gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_diag_flip_repr50.png"),res=450,units='in',width=7.75,height=4.5)
# pdf(paste0(imgpref, "epigenomes_gwas_flathg_top", SHOWGWAS, sub(".Rda", "",suffix), "_groupord_diag_flip_repr50.pdf"), width=7.75, height=4.5)
sp = 0.15
layout(matrix(c(1,2),2,1), heights=c(2,8), TRUE)
par(yaxs="i")
par(xaxs="i")
par(xpd=FALSE)
par(mar=c(0, 11, 3.75, 1.25))
meta.image(metamat[rnorder,], colvals=colvals, cex=0, horiz=F, useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
metaclass = sapply(colnames(metamat), capitalize)
text(y=seq(0,1, length.out=ncol(metamat)),
     x=par()$usr[2]+0.001*(par()$usr[2]-par()$usr[1]), 
     labels=metaclass, srt=0, adj=0, xpd=TRUE, cex=.3)
# Plot labels:
box.pad = 0.02
xw = lablist[[1]]
xx = space.1d(xw, box.pad=box.pad)
xx = space.1d(xx, box.pad=box.pad)
rx = c(0.05, 0.15, 0.4, 0.45)
x = par()$usr[4]+rx*(diff(par()$usr[3:4]))
text(x=xx, y=x[4], labels=lablist[[2]], col=lablist[[3]],
     srt=90, adj=0, xpd=TRUE, cex=.4)
par(xpd=TRUE)
segments(x0=xw, y0=x[2], x1=xx, y1=x[3], col=lablist[[3]])
segments(x0=xw, y0=x[1], x1=xw, y1=x[2], col=lablist[[3]])
par(xpd=FALSE)
par(mar=c(sp, 11, sp, 1.25))
image(epimat[rnorder, rev(diagord)], axes=FALSE,col=colred, zlim=c(zmin, zmax), useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
# add wordlet text:
# par(xpd=TRUE)
par(xpd=NA)
rx=c(0.005,0.01, 0.03, 0.034, 0.035)
x = par()$usr[1]-rx*(diff(par()$usr[1:2]))
xl = wordlets
xc = wcol
xx = seq(0 + 1e-3, 1.25 - 1e-3, length.out=length(xl))
xw = wind / (length(diagord) +1)
text(y=xx, x=x[5], labels=xl,
     srt=0, adj=1, xpd=NA, cex=termcex, col=wcol)
segments(x0=x[1], y0=xw, x1=x[2], y1=xw, col=wcol, lwd=.5)
segments(x0=x[2], y0=xw, x1=x[3], y1=xx, col=wcol, lwd=.5)
segments(x0=x[3], y0=xx, x1=x[4], y1=xx, col=wcol, lwd=.5)
par(xpd=FALSE)
# Add rectangles
rll = calc.breaks.rect(hcls=rev(vcut), vcls=acutnamed, colset)
rectdf = rll[[1]]
vccols = rll[[2]]
rord = order(rectdf$y1)
rectdf = rectdf[rord,]
vccols = vccols[rord]
hbreaks = calc.breaks.acut(acutnamed)
abline(h=unique(c(rectdf$y1, rectdf$y2)),lty='dotted', lwd=.5, col='darkgrey')
abline(v=hbreaks,lty='dotted', lwd=.5, col='darkgrey')
# Add grey - ubq rectangle:
rectdf = rbind(c(x1=par()$usr[1], x2=par()$usr[2],
                 y1=rectdf$y2[nrow(rectdf)], y2=par()$usr[4]), rectdf)
vccols = c('grey',vccols)
rect(xleft=rectdf$x1, xright=rectdf$x2,
     ybottom=rectdf$y1, ytop=rectdf$y2,
     border=vccols, lwd=.75)
par(xpd=NA)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
dev.off()

# -------------------------
# Count number enrichments:
# -------------------------
bmat = emat
bmat[bmat > 3] = 1
odf.nn = odf[which(odf$GROUP != 'None'),]
gr.mat = as.character(meta[rownames(emat), 'GROUP'])
tform = make.tform(gr.mat, u=odf.nn$GROUP)
bmat = t(bmat) %*% tform
count.df = data.frame(uid=rownames(bmat), 
                      total=apply(bmat, 1, sum), 
                      groups=apply(bmat > 0, 1, sum))
rownames(count.df) = count.df$uid
count.ord = names(sort(apply(bmat, 1, sum), decreasing=T))
count.df = count.df[count.ord,]
# Binary:
b.bmat = bmat[count.ord,] > 0
col.bmat = sweep(b.bmat, 2, 1:ncol(b.bmat), '*')

# Categorize:
norm.emat = sweep(emat, 2, apply(emat, 2, sum), '/')
toppct.emat = apply(norm.emat, 2, max)
barplot(toppct.emat[full.diagord])
which(full.vcut == 0)

# Multifact:
typedf = data.frame(uid=full.diagord, frac=toppct.emat[full.diagord], vcut=full.vcut, type=1)
rownames(typedf) = typedf$uid
typedf$type[typedf$frac < 0.5] = 2
typedf$type[typedf$vcut == 0] = 3
lgtype = paste0(c('Uni-factorial','Multi-factorial','Poly-factorial'), ' (', table(typedf$type),')')
coltype = c('indianred','steelblue','grey')

# Order by counts:
count.ord = names(sort(apply(bmat > 0, 1, sum), decreasing=F))
pdf(paste0(imgpref,'marginal_counts',sub(".Rda", ".pdf", suffix)), width=3, height=7)
layout(matrix(c(1:4), ncol=4), widths=c(.25, 1,1,1), TRUE)
sp = 0.15
rsp = 3
par(yaxs='i')
par(xaxs='i')
par(mar=c(rsp, 0, sp, sp))
image(t(as.matrix(typedf[count.ord, 'type'])), col=coltype, axes=F, ylab='', xlab='', useRaster=T)
par(mar=c(rsp, 0, sp, sp))
barplot(t(bmat[count.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab= '', yaxt='n', horiz=T)
legend('bottomright', lgtype, col=coltype, pch=15, cex=.6,pt.cex=1.25, inset=.01, box.col=NA, border='black')
par(mar=c(rsp, 0, sp, sp))
barplot(t(b.bmat[count.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab='', yaxt='n', horiz=T)
par(mar=c(rsp, 0, sp, sp))
image(t(col.bmat[count.ord,]), col=c('white', odf.nn$COLOR), axes=F, ylab='', xlab='',useRaster=T)
dev.off()

# Order as a heatmap:
gwas.ord = rev(full.diagord)
pdf(paste0(imgpref,'marginal_counts_bygroup',sub(".Rda", ".pdf", suffix)), width=3, height=7)
layout(matrix(c(1:4), ncol=4), widths=c(.25, 1,1,1), TRUE)
sp = 0.15
rsp = 3
par(yaxs='i')
par(xaxs='i')
par(mar=c(rsp, 0, sp, sp))
image(t(as.matrix(typedf[gwas.ord, 'type'])), col=coltype, axes=F, ylab='', xlab='', useRaster=T)
par(mar=c(rsp, 0, sp, sp))
barplot(t(bmat[gwas.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab= '', yaxt='n', horiz=T)
legend('bottomright', lgtype, col=coltype, pch=15, cex=.6,pt.cex=1.25, inset=.01, box.col=NA, border='black')
par(mar=c(rsp, 0, sp, sp))
barplot(t(b.bmat[gwas.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab='', yaxt='n', horiz=T)
par(mar=c(rsp, 0, sp, sp))
image(t(col.bmat[gwas.ord,]), col=c('white', odf.nn$COLOR), axes=F, ylab='', xlab='', useRaster=T)
dev.off()

pscount=0.5
bscale = log2(apply(bmat, 1, sum) + pscount)
l.bmat = sweep(bmat, 1, bscale / apply(bmat, 1, sum), '*')

# ONLY 1 tissue:
only1 = scan('fdrp1_names_only_1.tsv', 'c', sep='\n')
multi1 = scan('fdrp1_names_75over.tsv', 'c', sep='\n')
multi2 = scan('fdrp1_names_50over.tsv', 'c', sep='\n')
poly1 = scan('fdrp1_names_9overgroup.tsv', 'c', sep='\n')
poly2 = scan('fdrp1_names_17overgroup.tsv', 'c', sep='\n')
# Make matrix of these to show MK.

# Order as a heatmap:
gwas.ord = rev(full.diagord)
pdf(paste0(imgpref,'marginal_counts_bygroup_logscale',sub(".Rda", ".pdf", suffix)), width=3, height=7)
layout(matrix(c(1:4), ncol=4), widths=c(.25, 1,1,1), TRUE)
sp = 0.15
rsp = 3
par(yaxs='i')
par(xaxs='i')
par(mar=c(rsp, 0, sp, sp))
image(t(as.matrix(typedf[gwas.ord, 'type'])), col=coltype, axes=F, ylab='', xlab='', useRaster=T)
par(mar=c(rsp, 0, sp, sp))
barplot(t(l.bmat[gwas.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab= '', yaxt='n', horiz=T, xaxt='n')
lbpar= c(0,2^(0:9))
atpar = log2(lbpar + pscount)
axis(1, at=atpar, labels=lbpar,las=2, cex.axis=.5, lwd=.5)
abline(v=atpar, lwd=.25, col='grey50', lty='dashed')
barplot(t(l.bmat[gwas.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab= '', yaxt='n', horiz=T, xaxt='n', add=TRUE)
par(mar=c(rsp, 0, sp, sp))
barplot(t(b.bmat[gwas.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab='', yaxt='n', horiz=T)
legend('bottomright', lgtype, col=coltype, pch=15, cex=.6,pt.cex=1.25, inset=.01, box.col=NA, border='black')
par(mar=c(rsp, 0, sp, sp))
image(t(col.bmat[gwas.ord,]), col=c('white', odf.nn$COLOR), axes=F, ylab='', xlab='', useRaster=T)
dev.off()



# ----------------------------------
# Count number enrichments FOR TREE:
# ----------------------------------
eprefix = 'enhancers_e2500_'
apref = 'cons_parent'
snpfile = paste0(regdir, eprefix, apref, '_logreg_all_wsnp_adj1000_1.Rda')
load(snpfile, verbose=T)
# Loads: rmat smat nlist
load(paste0(eprefix, 'kept_allgwas_ordered_adj1000_1.Rda'))
print(length(Znam))  # Traits to look at/plot
kuid = rownames(rmat)
kuid = kuid[kuid%in% Znam]

bmat = t(rmat[kuid,])
bmat[bmat > 3] = 1
nodetissue = nodetissue[order(nodetissue$node),]
gr.mat = nodetissue$GROUP
tform = make.tform(gr.mat, u=odf$GROUP)
bmat = t(bmat) %*% tform
count.df = data.frame(uid=rownames(bmat), 
                      total=apply(bmat, 1, sum), 
                      groups=apply(bmat > 0, 1, sum))
rownames(count.df) = count.df$uid
count.ord = names(sort(apply(bmat, 1, sum), decreasing=T))
count.df = count.df[count.ord,]
# Binary:
b.bmat = bmat[count.ord,] > 0
col.bmat = sweep(b.bmat, 2, 1:ncol(b.bmat), '*')

# Order by counts:
count.ord = names(sort(apply(bmat > 0, 1, sum), decreasing=F))
pdf(paste0(imgpref,'tree_marginal_counts',sub(".Rda", ".pdf", suffix)), width=3, height=7)
layout(matrix(c(1:3), ncol=3), TRUE)
sp = 0.15
rsp = 3
par(yaxs='i')
par(xaxs='i')
par(mar=c(rsp, 0, sp, sp))
barplot(t(bmat[count.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab= '', yaxt='n', horiz=T)
par(mar=c(rsp, 0, sp, sp))
barplot(t(b.bmat[count.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab='', yaxt='n', horiz=T)
par(mar=c(rsp, 0, sp, sp))
image(t(col.bmat[count.ord,]), col=c('white', odf.nn$COLOR), axes=F, ylab='', xlab='',useRaster=T)
dev.off()

# Order as a heatmap:
gwas.ord = rev(Znam)[rev(Znam) %in% kuid]
pdf(paste0(imgpref,'tree_marginal_counts_bygroup',sub(".Rda", ".pdf", suffix)), width=3, height=7)
layout(matrix(c(1:3), ncol=3), TRUE)
sp = 0.15
rsp = 3
par(yaxs='i')
par(xaxs='i')
par(mar=c(rsp, 0, sp, sp))
barplot(t(bmat[gwas.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab= '', yaxt='n', horiz=T)
par(mar=c(rsp, 0, sp, sp))
barplot(t(b.bmat[gwas.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab='', yaxt='n', horiz=T)
par(mar=c(rsp, 0, sp, sp))
image(t(col.bmat[gwas.ord,]), col=c('white', odf.nn$COLOR), axes=F, ylab='', xlab='', useRaster=T)
dev.off()

pscount=0.5
bscale = log2(apply(bmat, 1, sum) + pscount)
l.bmat = sweep(bmat, 1, bscale / apply(bmat, 1, sum), '*')

# Order as a heatmap:
gwas.ord = rev(Znam)[rev(Znam) %in% kuid]
pdf(paste0(imgpref,'tree_marginal_counts_bygroup_logscale',sub(".Rda", ".pdf", suffix)), width=3, height=7)
layout(matrix(c(1:3), ncol=3), TRUE)
sp = 0.15
rsp = 3
par(yaxs='i')
par(xaxs='i')
par(mar=c(rsp, 0, sp, sp))
barplot(t(l.bmat[gwas.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab= '', yaxt='n', horiz=T, xaxt='n')
lbpar= c(0,2^(0:9))
atpar = log2(lbpar + pscount)
axis(1, at=atpar, labels=lbpar,las=2, cex.axis=.5, lwd=.5)
abline(v=atpar, lwd=.25, col='grey50', lty='dashed')
barplot(t(l.bmat[gwas.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab= '', yaxt='n', horiz=T, xaxt='n', add=TRUE)
par(mar=c(rsp, 0, sp, sp))
barplot(t(b.bmat[gwas.ord,]), col=odf.nn$COLOR, border=NA, ylab='', xlab='', yaxt='n', horiz=T)
par(mar=c(rsp, 0, sp, sp))
image(t(col.bmat[gwas.ord,]), col=c('white', odf.nn$COLOR), axes=F, ylab='', xlab='', useRaster=T)
dev.off()

sum(kuid %in% full.diagord)






# ------------------------------------------------
# Measure how many enrichments vs. old epigenomes:
# ------------------------------------------------
gwdf = data.frame(t(emat))
gwdf$uid = rownames(gwdf)
gwdf = gather(gwdf, id, logpval, -uid)
gwdf = gwdf[gwdf$logpval > 0,]
# Mapping:
passdf = merge(gwdf, meta[,c('id','GROUP','infoline','Project')])

oldepi = c('ENCODE 2012','Roadmap 2015')
oldpassdf = filter(passdf, Project %in% oldepi)
# Raw # of enrichments + %
print(paste(nrow(oldpassdf),"to",nrow(passdf)))
# New gwas studies (PMIDxGWAS) annotated # + %
oldst = unique(oldpassdf$uid)
newst = unique(c(passdf$uid, oldpassdf$uid))
print(paste(length(oldst), "to", length(newst)))
idmap = meta[,c('id','Project', 'GROUP','infoline')]

topenr = c()
topnum = c()
for (st in newst){
    stdf = filter(passdf, uid == st)
    mrow = stdf[which.max(stdf$logpval),]
    sid = as.character(mrow$id)
    topenr = c(topenr, sid)
    topnum = c(topnum, mrow$logpval)
}
topdf = data.frame(id=topenr, num=topnum, uid=newst)
topdf = merge(topdf, idmap)

tdf = aggregate(uid ~ Project, topdf, length)
tdf = tdf[order(tdf$uid, decreasing=T),]
tcols = colvals[['project']][tdf$Project]
vals = tdf$uid
names(vals) = tdf$Project
tdf$lab = paste0(tdf$Project, ' (', tdf$uid, ')')

png(paste0(imgpref,'numtop_byproj_cutoff',sub(".Rda", ".png", suffix)),res=450,units='in',width=6,height=2)
# pdf(paste0(imgpref,'numtop_byproj_cutoff',sub(".Rda", ".pdf", suffix)), width=6,height=2)
par(mar=c(3,.5,1,3.75))
par(yaxs='i')
par(xaxs='i')
barplot(tdf$uid, width=.5, col=tcols, # names.arg=tdf$Project, 
        xlab='', space=0.1, xlim = c(0, (round(max(tdf$uid) / 50) + 1) * 50),
        main='', border=NA, horiz=T)
abline(v=seq(0,max(tdf$uid), 50), col='white', lwd=1, lty='dotted')
abline(v=seq(0,max(tdf$uid), 100), col='white', lwd=2)
mtext('# Strongest enrichment', side=3, line=0, cex=1.2)
mtext('Number of GWAS', side=1, line=2, cex=1)
txtoffset = diff(par()$usr[3:4]) / (2 * nrow(tdf))
txtloc = seq(par()$usr[3] + txtoffset, par()$usr[4] - txtoffset, length.out=nrow(tdf))
text(y=txtloc, x=tdf$uid + 3, tdf$lab, 
     xpd=TRUE, srt=0, adj=0, cex=1)
dev.off()


# -------------------------------
# Also count gwas by new vs. not:
# -------------------------------
epidf = rbind(data.frame(Project=c('ENCODE 2012','Roadmap 2015'), pstatus='Old'),
              data.frame(Project=c('ENCODE (New)','Roadmap (New)', 'GGR'), pstatus='New'))
tdf = aggregate(uid ~ pstatus + Project, merge(topdf, epidf), length)
tdf = tdf[order(tdf$uid, decreasing=F),]
vals = tdf$uid
names(vals) = tdf$pstatus
tdf$lab = paste0(tdf$pstatus, ' (', tdf$uid, ')')

twide = spread(tdf[,c('pstatus','Project','uid')], Project, uid)
tmat = as.matrix(twide[,-1])
tmat[is.na(tmat)] = 0
rownames(tmat) = twide$pstatus
ylim = c(0, (round(max(apply(tmat, 1, sum)) / 50) + 1) * 50)

png(paste0(imgpref,'numtop_byprojstatus_cutoff',sub(".Rda", ".png", suffix)),res=450,units='in',width=2,height=6)
# pdf(paste0(imgpref,'numtop_byprojstatus_cutoff',sub(".Rda", ".pdf", suffix)), width=2,height=6)
par(mar=c(2,3,1,.5))
par(yaxs='i')
par(xaxs='i')
barplot(t(tmat), width=.5, col=tcols[colnames(tmat)],
        xlab='', space=0.1, ylim=ylim, main='', border=NA)
# abline(h=seq(50,max(tdf$uid), 50), col='white', lwd=1, lty='dotted')
abline(h=seq(100,max(tdf$uid), 100), col='white', lwd=1)
mtext('# Strongest\nenrichment', side=3, line=-1.5, cex=1.2)
mtext('Number of GWAS', side=2, line=2, cex=1)
txtoffset = diff(par()$usr[1:2]) / (2 * nrow(tdf))
txtloc = seq(par()$usr[1] + txtoffset, par()$usr[2]-txtoffset, length.out=nrow(tdf))
text(x=txtloc, y=tdf$uid + 10, tdf$lab, xpd=TRUE, srt=0, adj=.5, cex=.85)
dev.off()



# UpSet plot (pseudo venn):
redpassdf = unique(passdf[,c('uid','Project')])
ulist = list()
for (pset in unique(redpassdf$Project)){
    ulist[[pset]] = unique(redpassdf$uid[redpassdf$Project == pset])
}

png(paste0(imgpref,'upset_byproj',sub(".Rda", ".png", suffix)),res=450,units='in',width=9,height=6)
upset(fromList(ulist), order.by = "freq")
dev.off()

pdf(paste0(imgpref,'upset_byproj',sub(".Rda", ".pdf", suffix)),width=9,height=6)
upset(fromList(ulist), order.by = "freq")
dev.off()


m2 = make_comb_mat(ulist)
orig.ord = names(set_size(m2))
# ord = names(sort(set_size(m2)))
ord = rev(c("ENCODE (New)", "Roadmap (New)", "Roadmap 2015", "ENCODE 2012", "GGR"))
scol = colvals[['project']][orig.ord]
names(scol) = NULL

# png(paste0(imgpref,'upset2_byproj',sub(".Rda", ".png", suffix)),res=450,units='in',width=9,height=3)
pdf(paste0(imgpref,'upset2_byproj',sub(".Rda", ".pdf", suffix)), width=9,height=3)
UpSet.v2(m2, set_col=scol, set_order=ord,
         comb_order=order(-comb_size(m2)), bg_pt_col="#F0F0F0")
dev.off()

# Rarefaction-style curves (both types)
rdf = passdf
satord = c()
satnum = c()
sattr = c()
while (nrow(rdf) != 0){
    ag = aggregate(uid ~ id, rdf, length)
    mrow = ag[which.max(ag$uid),]
    sid = as.character(mrow[[1]])
    npm = mrow[[2]]
    satord = c(satord, sid)
    satnum = c(satnum, npm)
    # Find top enrichment:
    ind = which(rdf$id == sid)
    topind = ind[order(rdf$logpval[ind], decreasing=T)]
    toptr = sub("^[0-9]* - ", "", rdf$uid[head(topind, 6)])
    sattr = c(sattr, paste(toptr, collapse=', '))
    # Remove:
    pid = filter(rdf, id == sid)$uid
    rdf = filter(rdf, !(uid %in% pid))
}
satdf = data.frame(id=satord, num=satnum, cs=cumsum(satnum), trait=sattr)
satdf = merge(satdf, idmap)
satdf = satdf[order(satdf$cs),]
satdf$label = paste0(satdf$infoline, " (", satdf$trait, ")")

satcols = colvals[['project']][satdf$Project]
satgroupcols = colvals[['group']][satdf$GROUP]
# png(paste0(imgpref,'saturation_all_cutoff', sub(".Rda", ".png", suffix)),res=450,units='in',width=7,height=4)
pdf(paste0(imgpref,'saturation_all_cutoff', sub(".Rda", ".pdf", suffix)),width=7,height=4)
par(mar=c(3,4.5,.15,.15))
plot(1:nrow(satdf), satdf$cs, pch=19, col=satcols, ylim=c(0, max(satdf$cs)), 
     cex=.5, axes=F, ylab='', xlab='')
mtext('Number of epigenomes', 1, cex=.85, line=1.75)
mtext('Cumulative number\nof annotated GWAS', 2, cex=.85, line=2.5)
abline(h=seq(0,max(satdf$cs),100), col='darkgrey',lty='dotted')
box(lwd=.5)
nlabel = 25
cex = 0.5
labels = satdf$infoline[1:nlabel]
# labels = satdf$label[1:nlabel]
textHeight <- graphics::strheight(labels, cex = cex)
textWidth <- graphics::strwidth(labels, cex = cex)
box.pad = textHeight[1] * 1.15
y1 = (satdf$cs[1:nlabel] - 0.005 * diff(par()$usr[3:4]))
x1 = 1:nlabel + 0.0025 * diff(par()$usr[1:2])
x2 = 1:nlabel + ((1:nlabel > 10) * 0.02 +  0.04) * diff(par()$usr[1:2])
x3 = 1:nlabel + ((1:nlabel > 10) * 0.02 +  0.0425) * diff(par()$usr[1:2])
yw = seq(satdf$cs[1]-75, satdf$cs[nlabel]+5,length.out=nlabel)
yy = space.1d(yw, box.pad=box.pad, lim=c(-15, max(satdf$cs)))
yy = space.1d(yy, box.pad=box.pad, lim=c(-15, max(satdf$cs)))
# Add box:
rect(xleft=x3 - .1*textWidth/2, xright=x3 + 2.1*textWidth/2, 
     ybottom=yy - 1.1*textHeight/2, ytop=yy + 1.1*textHeight/2,
     col="white", border=NA)
segments(x0=x1, y0=y1, x1=x2, y1=yy, col=satgroupcols[1:nlabel], lwd=.5)
segments(x0=x2, y0=yy, x1=x3, y1=yy, col=satgroupcols[1:nlabel], lwd=.5)
text(x=x3, y=yy,  labels=labels, 
     srt=0, adj=0, cex=cex, col=satgroupcols[1:nlabel])
legend('bottomright', names(colvals[['project']]), col=colvals[['project']], 
       pch=19, cex=.75, inset=.01, box.col=NA, title=expression(bold('Project')))
cex.axis=.75
axis(1, lwd=.5, las=1, cex.axis=cex.axis, padj=-1)
axis(2, lwd=.5, las=1, cex.axis=cex.axis, hadj=1)
dev.off()


satdf$trait.label = sapply(satdf$trait, width = 40, split.text)
satcols = colvals[['project']][satdf$Project]
satgroupcols = colvals[['group']][satdf$GROUP]
png(paste0(imgpref,'saturation_all_wlab_cutoff', sub(".Rda", ".png", suffix)),res=450,units='in',width=7,height=4)
# pdf(paste0(imgpref,'saturation_all_wlab_cutoff', sub(".Rda", ".pdf", suffix)),width=7,height=4)
par(mar=c(3,4.5,.15,.15))
plot(1:nrow(satdf), satdf$cs, pch=19, col=satcols, ylim=c(0, max(satdf$cs)), 
     cex=.5, axes=F, ylab='', xlab='')
mtext('Number of epigenomes', 1, cex=.85, line=1.75)
mtext('Cumulative number\nof annotated GWAS', 2, cex=.85, line=2.5)
abline(h=seq(0,max(satdf$cs),100), col='darkgrey',lty='dotted', cex=.5)
box(lwd=.5)
nlabel = 25
cex = 0.5
labels = satdf$infoline[1:nlabel]
# labels = satdf$label[1:nlabel]
lab.textHeight <- graphics::strheight(labels, cex = cex)
lab.textWidth <- graphics::strwidth(labels, cex = cex)
box.pad = lab.textHeight[1] * 1.15
y1 = (satdf$cs[1:nlabel] - 0.005 * diff(par()$usr[3:4]))
x1 = 1:nlabel + 0.0025 * diff(par()$usr[1:2])
x2 = 1:nlabel + ((1:nlabel > 10) * 0.02 +  0.04) * diff(par()$usr[1:2])
x3 = 1:nlabel + ((1:nlabel > 10) * 0.02 +  0.0425) * diff(par()$usr[1:2])
yw = seq(satdf$cs[1]-75, satdf$cs[nlabel]+5,length.out=nlabel)
yy = space.1d(yw, box.pad=box.pad, lim=c(-15, max(satdf$cs)))
yy = space.1d(yy, box.pad=box.pad, lim=c(-15, max(satdf$cs)))
ntrait = 10
tr.cex = 0.4
traits = satdf$trait.label[1:ntrait]
textHeight <- graphics::strheight(traits, cex = tr.cex)
textWidth <- graphics::strwidth(traits, cex = tr.cex)
# Trait segments:
tx1 = (x3 + lab.textWidth)[1:ntrait]
tx2 = par()$usr[1] + 0.45 * diff(par()$usr[1:2])
tx3 = par()$usr[1] + (0.77 - 0.24 * (1:ntrait > ntrait/2)) * diff(par()$usr[1:2])
tx4 = par()$usr[1] + (0.78 - 0.24 * (1:ntrait > ntrait /2 )) * diff(par()$usr[1:2])
ty1 = yy[1:ntrait]
tyy = c(seq(max(satdf$cs) * 0.05, max(satdf$cs) * .55,length.out=ntrait / 2), 
        seq(max(satdf$cs) * 0.3, max(satdf$cs) * .83,length.out=ntrait / 2))
segments(x0=tx1, y0=ty1, x1=tx2, y1=ty1, col=satgroupcols[1:ntrait], lwd=.5)
segments(x0=tx2, y0=ty1, x1=tx3, y1=tyy, col=satgroupcols[1:ntrait], lwd=.5)
rect(xleft=tx4 - .02*textWidth/2, xright=tx4 + 2*textWidth/2, 
     ybottom=tyy - 1.1*textHeight/2, ytop=tyy + 1.1*textHeight/2,
     col="white", border=NA)
segments(x0=tx3, y0=tyy, x1=tx4, y1=tyy, col=satgroupcols[1:ntrait], lwd=.5)
text(x=tx4, y=tyy,  labels=traits, srt=0, adj=0, cex=tr.cex, col=satgroupcols[1:ntrait])
# Add the labels last (goes over lines)
rect(xleft=x3 - .02*lab.textWidth/2, xright=x3 + 2.02*lab.textWidth/2, 
     ybottom=yy - 1.1*lab.textHeight/2, ytop=yy + 1.1*lab.textHeight/2,
     col="white", border=NA)
segments(x0=x1, y0=y1, x1=x2, y1=yy, col=satgroupcols[1:nlabel], lwd=.5)
segments(x0=x2, y0=yy, x1=x3, y1=yy, col=satgroupcols[1:nlabel], lwd=.5)
text(x=x3, y=yy,  labels=labels, srt=0, adj=0, cex=cex, col=satgroupcols[1:nlabel])
legend('topleft', names(colvals[['project']]), col=colvals[['project']], 
       pch=19, cex=.5, inset=.01, box.col=NA, title=expression(bold('Project')))
cex.axis=.75
axis(1, lwd=.5, las=1, cex.axis=cex.axis, padj=-1)
axis(2, lwd=.5, las=1, cex.axis=cex.axis, hadj=1)
dev.off()



rdf = passdf
mdf = aggregate(logpval ~ uid, rdf, max)
rdf = merge(mdf, rdf)
satord = c()
satnum = c()
while (nrow(rdf) != 0){
    ag = aggregate(uid ~ id, rdf, length)
    mrow = ag[which.max(ag$uid),]
    sid = as.character(mrow[[1]])
    npm = mrow[[2]]
    satord = c(satord, sid)
    satnum = c(satnum, npm)
    # Remove:
    pid = filter(rdf, id == sid)$uid
    rdf = filter(rdf, !(uid %in% pid))
}
satdf = data.frame(id=satord, num=satnum, cs=cumsum(satnum))
satdf = merge(satdf, idmap)
satdf = satdf[order(satdf$cs),]

satcols = colvals[['project']][satdf$Project]
satgroupcols = colvals[['group']][satdf$GROUP]

# png(paste0(imgpref,'saturation_all_max_cutoff', sub(".Rda", ".png", suffix)),res=450,units='in',width=7,height=4)
pdf(paste0(imgpref,'saturation_all_max_cutoff', sub(".Rda", ".pdf", suffix)),width=7,height=4)
par(mar=c(3,4.5,.15,.15))
plot(1:nrow(satdf), satdf$cs, pch=19, col=satcols, ylim=c(0, max(satdf$cs)), 
     cex=.3, axes=F, ylab='', xlab='')
mtext('Number of epigenomes', 1, cex=.85, line=1.75)
mtext('Cumulative number\nof annotated GWAS', 2, cex=.85, line=2.5)
abline(h=seq(0,max(satdf$cs),100), col='darkgrey',lty='dotted')
box(lwd=.5)
nlabel = 50
# lines:
cex = 0.25
textHeight <- graphics::strheight(satdf$infoline[1:nlabel], cex = cex)
textWidth <- graphics::strwidth(satdf$infoline[1:nlabel], cex = cex)
box.pad = textHeight[1] * 1.15
y1 = (satdf$cs[1:nlabel] - 0.005 * diff(par()$usr[3:4]))
x1 = 1:nlabel + 0.0025 * diff(par()$usr[1:2])
x2 = 1:nlabel + ((1:nlabel > 25) * 0.02 +  0.04) * diff(par()$usr[1:2])
x3 = 1:nlabel + ((1:nlabel > 25) * 0.02 +  0.05) * diff(par()$usr[1:2])
yw = seq(satdf$cs[1]-15, satdf$cs[nlabel]+25,length.out=nlabel)
yy = space.1d(yw, box.pad=box.pad, lim=c(-15, max(satdf$cs)))
yy = space.1d(yy, box.pad=box.pad, lim=c(-15, max(satdf$cs)))
# Add box:
rect(xleft=x3 - .1*textWidth/2, xright=x3 + 2.1*textWidth/2, 
     ybottom=yy - 1.1*textHeight/2, ytop=yy + 1.1*textHeight/2,
     col="white", border=NA)
segments(x0=x1, y0=y1, x1=x2, y1=yy, col=satgroupcols[1:nlabel], lwd=.5)
segments(x0=x2, y0=yy, x1=x3, y1=yy, col=satgroupcols[1:nlabel], lwd=.5)
text(x=x3, y=yy,  labels=satdf$infoline[1:nlabel], 
     srt=0, adj=0, cex=cex, col=satgroupcols[1:nlabel])
legend('bottomright', names(colvals[['project']]), col=colvals[['project']], 
       pch=19, cex=.75, inset=.01, box.col=NA, title=expression(bold('Project')))
cex.axis=.75
axis(1, lwd=.5, las=1, cex.axis=cex.axis, padj=-1)
axis(2, lwd=.5, las=1, cex.axis=cex.axis, hadj=1)
dev.off()





# Rarefaction by project:
projord = c("ENCODE 2012", "Roadmap 2015", 
            "Roadmap (New)", "ENCODE (New)", "GGR")
# Rarefaction-style curves (both types)
rdf = passdf
satord = c()
satnum = c()
topproj = projord[1]
while (nrow(rdf) != 0){
    subdf = filter(rdf, Project == topproj)
    if (nrow(subdf) == 0){
        projord = projord[-1]
        topproj = projord[1]
        subdf = filter(rdf, Project == topproj)
    }
    ag = aggregate(uid ~ id, subdf, length)
    mrow = ag[which.max(ag$uid),]
    sid = as.character(mrow[[1]])
    npm = mrow[[2]]
    satord = c(satord, sid)
    satnum = c(satnum, npm)
    # Remove:
    pid = filter(rdf, id == sid)$uid
    rdf = filter(rdf, !(uid %in% pid))
}
satdf = data.frame(id=satord, num=satnum, cs=cumsum(satnum))
satdf = merge(satdf, idmap)
satdf = satdf[order(satdf$cs),]

# Which NEW TISSUES are the most informative?
satcols = colvals[['project']][satdf$Project]
satgroupcols = colvals[['group']][satdf$GROUP]
png(paste0(imgpref,'saturation_byproj_cutoff',sub(".Rda", ".png", suffix)),res=300,units='in',width=7,height=12)
par(mar=c(4.5,4.5,1,1))
plot(satdf$cs, rev(1:nrow(satdf)), pch=19, cex=.5, col=satcols, xlim=c(0, max(satdf$cs)), 
     ylab='Number of epigenomes', xlab='Cumulative number of GWAS')
abline(v=seq(0,max(satdf$cs),50), col='darkgrey',lty='dotted')
text(satdf$cs - 0.015 * diff(par()$usr[1:2]), rev(1:nrow(satdf)), labels=satdf$infoline, 
     srt=0, adj=1, cex=.4, col=satgroupcols)
legend('bottom', names(colvals[['project']]), col=colvals[['project']], 
       pch=19,inset=.01, box.col=NA, title=expression(bold('Project')))
legend('bottomleft', names(colvals[['group']]), col=colvals[['group']], 
       pch=19,inset=.01, box.col=NA, title=expression(paste(bold("Group"))))
dev.off()




