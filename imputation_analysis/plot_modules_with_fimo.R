#!/usr/bin/R
# -----------------------------------------------------
# Plot modules from CHMM/other clustering with fimo
# -----------------------------------------------------
# TODO: NOTE: THE BACKGROUND FOR EPIGENOMES IS NOT APPROPRIATE. 
# FIXME: DOUBLECOUNTING - ESPECIALLY PROBLEMATIC IN EMBRYONIC REGIONS
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
fset = 'q10'
fimofile = paste0(fset, '_sig_merged_fimo_pval.RData')
filepref = 'cls_merge2_wH3K27ac100_300'
tagline = 'ChromHMM Enhancers'
# Or directory?
# filepref = 'cls_merge2_wH3K27ac100_raw'
# tagline = 'ChromHMM Enhancers (on Epigenomes)'
# filepref = 'prom_cls/cls_merge2_wH3K27ac100_300'
# tagline = 'ChromHMM Promoters'
imgdir = paste0(img, "clusters/") 

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need clusters filename")
} else {        
    fimofile = args[1]
    filepref = args[2]
    tagline = args[3]
    fset = args[4]
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


# Prefix / plot directory:
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'clusters_fimo_', sub("/","_", filepref), '_', fset, '_')
print(paste("[STATUS] Plotting under prefix", imgpref))

# ---------------------------------
# Load all of the fimo data (RDATA)
# ---------------------------------
load(fimofile)
print(head(sigdf))
sigdf$cls = paste0('c', sigdf$id)
ids = sort(unique(sigdf$cls))

# Filter for log2fc > 0.5:
fccutoff = 0.5
sigdf = filter(sigdf, abs(log2fc) >= fccutoff)

# Filter bad motifs:
motifs = sort(unique(sigdf$motif))
outsp = c('Mus_musculus', 'Tetraodon_nigroviridis',
          'Drosophila_melanogaster','Xenopus_laevis')
remove.mot = unlist(sapply(outsp, function(x){grep(x, motifs)}))
keep.motifs = as.character(motifs[-remove.mot])
sigdf = filter(sigdf, motif %in% keep.motifs)

# Reduce motif names:
motdf = unique(sigdf[,c('motif','db')])
patternlist=c("_MA.*JASPAR", "_M[0-9].*_1.02_CIS-BP", "_HUMAN.H11MO.*_HUMAN", 
              "_full.*_EUKARYOTE", "_DBD.*_EUKARYOTE", "^_")
redmot = as.character(motdf$motif)
for (pat in patternlist){ redmot = sub(pat, "",redmot) }
motdf$rmotif = toupper(redmot)
motdf$db = as.character(motdf$db)
motdf$dmotif = paste0(motdf$rmotif, " (", motdf$db, ")")
# add to sigdf:
sigdf = merge(sigdf, motdf)

# TODO: Merge motif families

# Wide matrix:
siglong = aggregate(log2fc ~ cls + rmotif, sigdf, max)
wide = spread(siglong, cls, log2fc, fill=0)
smat = as.matrix(wide[,-1])
rownames(smat) = wide$rmotif

# TODO: Q: Keep motifs with at least log2 fc = 2?
# NOTE: Maybe for the cls, not for epi - not enough signal - diluted.
smarg = apply(smat, 1,function(x){max(abs(x))})
# keep.marg = names(which(smarg > 2))
# smat = smat[keep.marg,]

# Check if epigenomes:
namesfile = paste0(filepref, '_names.tsv')
epinames = scan(namesfile,'c')
is.epi = (length(grep("BSS",epinames)) == ncol(smat) || 
          sub(".*_","", filepref) == 'raw')
if (is.epi){
    print("Is enrichment on epigenomes, plotting with epigenomes")
    clsn = paste0('c',1:length(epinames) - 1)
    names(epinames) = clsn
    colnames(smat) = epinames[colnames(smat)]
}

# ----------------------------
mins = min(smat[!is.infinite(smat)])
maxs = max(smat[!is.infinite(smat)])
# Raw values, not thresholded:
holdmat = smat
holdmat[holdmat < mins] = mins
holdmat[holdmat > maxs] = maxs

# Thresholded 
tcut = 2
smat[smat > tcut] = tcut
smat[smat < -tcut] = -tcut

# Reorder alone:
method='euclidean'
r2 = reord(t(holdmat),method)
r3 = reord(holdmat, method)
epiord = rownames(r2)
motord = rownames(r3)
finalmat = t(smat[motord, epiord])
keep.cls = list()
keep.cls[[tagline]] = rownames(finalmat)


plot.fimo <-  function(fmat, set, ordll, tcut=2, ablwd=1,
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
    ll = diag.mat(tmat, ratio=0.25)
    tmp = ll[[1]]
    cto = ll[[3]] # For breaks
    diagord = colnames(tmp)
    gmat = gmat[,diagord]
    vcut =  c(cto[cto ==0], clscuts[cto])
    hbreaks = calc.breaks.acut(vcut)
    # Threshold for plotting:
    if (!calc.only){
        # TODO: check not flipped?
        image(gmat, axes=FALSE,col=rev(colrb), zlim=c(-tcut, tcut))
        abline(h=par()$usr[3:4],lty=1,lw=0.5)
        abline(v=par()$usr[1:2],lty=1,lw=0.5)
        abline(v=vbreaks,lty='dotted',lw=ablwd, col='darkgrey')
        abline(h=hbreaks,lty='dotted',lw=ablwd, col='darkgrey')
    }
    return(list(diagord, vcut))
}


# Set groupdist parameters
rnorder = rngroup # ord.gwas
acutnamed = acutgroup.nam
set = tagline 
lablist = lablist.group

dist.order = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, subset=TRUE, calc.only=TRUE)
gll = plot.fimo(finalmat, set, dist.order, motord=motord, calc.only=TRUE)
motnam = gll[[1]]

# Make plot:
png(paste0(imgpref,'groupord.png'),res=450,units='in',width=12,height=17)
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
image(finalmat[,motnam], col='white', axes=F)
text(y=seq(0,1, length.out=length(motnam)),
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=motnam, srt=0, adj=1, xpd=TRUE,cex=.25)
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
# FIMO plot
par(mar=c(0.25, 0.25, 0, 0.25))
gll = plot.fimo(finalmat, set, dist.order, motord=motord)
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
     border=vccols, lwd=1)
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
image(avail, axes=F, col=c('white', 'darkgrey'))
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
par(mar=c(.25, 0, 0, 0.25))
image(avail, axes=F, col='white')
text(x=seq(0,1, length.out=length(main.marks)),
     y=par()$usr[4], 
     labels=main.marks, srt=90, adj=1, xpd=TRUE,cex=.5)
dev.off()


# ---------------------
# Plot only the motifs: - 
# TODO: All cls (so that it matches other figures.
# ---------------------
# TODO: Must plot all cls!
png(paste0(imgpref,'groupord_onlyfimo.png'),res=450,units='in',width=7.75,height=15)
layout(matrix(c(1:2),1,2), widths=c(2,8), TRUE)
par(xpd=FALSE)
par(yaxs="i")
par(xaxs="i")
dist.order = list()
sp=0.15
# Motif labels:
par(mar=c(sp, 6, sp, sp))
par(xpd=NA)
image(finalmat[,motnam], col='white', axes=F)
text(y=seq(0,1, length.out=length(motnam)),
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=motnam, srt=0, adj=1, xpd=TRUE,cex=.25)
par(xpd=FALSE)
# Calculate centers:
set = tagline
dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, title=FALSE, subset=TRUE, calc.only=TRUE)
par(mar=c(sp, sp, sp, sp))
gll = plot.fimo(finalmat, set, dist.order[[set]], motord=motord)
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
     border=vccols, lwd=1)
par(xpd=NA)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[3], 
     ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
par(xpd=FALSE)
dev.off()


# Reduce to specific motifs:
finalmat = t(smat[motord, epiord])
kmot = names(which(apply(abs(finalmat) > .25, 2, sum) < nrow(finalmat) * .25))
kmot = kmot[nchar(kmot) < 15]
finalmat = finalmat[,kmot]
motord = colnames(finalmat)

dist.order = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, subset=TRUE, calc.only=TRUE)
gll = plot.fimo(finalmat, set, dist.order, motord=motord, calc.only=TRUE)
motnam = gll[[1]]
gcut = gll[[2]]
gcol = c('grey', colset)[gcut + 1]

png(paste0(imgpref,'groupord_onlyfimo_specific.png'),res=450,units='in',width=7.75,height=5)
layout(matrix(c(1:2),1,2), widths=c(1.7,8), TRUE)
par(xpd=FALSE)
par(yaxs="i")
par(xaxs="i")
dist.order = list()
sp = 0.15
# Motif labels:
par(mar=c(sp, 6, sp, sp))
par(xpd=NA)
image(finalmat[,motnam], col='white', axes=F)
# text(y=seq(0,1, length.out=length(motnam)),
#      x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
#      labels=motnam, srt=0, adj=1, xpd=TRUE,cex=.35)
nb = 2
for (i in c(1:nb - 1)){
    fm = 1:length(motnam)
    fm = floor(fm / 8)
    fm = which((fm + 1) %% nb == i )
    box.pad = 0.006
    xold = seq(0,1, length.out=length(motnam))[fm]
    xx = xold * 2
    # while (sum(xx != xold) > 0){
    #     xold = xx
        xx = space.1d(xx, box.pad=box.pad)
        xx = space.1d(xx, box.pad=box.pad)
    # }
    # xx = space.1d(xx, box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
    text(y=xx, x=par()$usr[2]- (0.001 + i * 3.5)*(par()$usr[2]-par()$usr[1]), 
         labels=motnam[fm], srt=0, adj=1, xpd=TRUE,cex=.2, col=gcol[fm])
}
par(xpd=FALSE)
# Calculate centers:
set = tagline
dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, title=FALSE, subset=TRUE, calc.only=TRUE)
rll = calc.breaks.rect(hcls=acutnamed, vcls=dist.order[[set]][[3]], colset)
clsrectdf = rll[[1]]
clsvccols = rll[[2]]
clsrectdf = rbind(c(x1=par()$usr[1], x2=clsrectdf$x1[1],
                 y1=par()$usr[3], y2=par()$usr[4]), clsrectdf)
clsvccols = c('grey',clsvccols)
par(mar=c(sp, 0, sp, sp))
gll = plot.fimo(finalmat, set, dist.order[[set]], motord=motord, ablwd=0.5)
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

