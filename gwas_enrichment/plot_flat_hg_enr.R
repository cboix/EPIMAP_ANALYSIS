#!/usr/bin/R
# ----------------------------------------------
# Plot the results of flat enrichments analysis:
# ----------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'auxiliary_gwastree_functions.R'))

library(GenomicRanges)
library(dplyr)
library(cba)
library(ggplot2)

# bedfile='/broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust/cls_merge2_wH3K27ac100_nonovl_300_seed1_assignments.fixed.bed'
# filepref='/broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust/cls_merge2_wH3K27ac100_nonovl_300_seed1'
# Arguments:
args=commandArgs(TRUE)
tagline=''
if (length(args)==0) {
    print("Using default arguments. Requires bedfile at least...")
} else {        
    bedfile = args[1]
    filepref = args[2]
    if (length(args) > 2){
        resultdir = args[3]
    } else {
        resultdir = sub("\\..*", "/", bedfile)
    }
    if (length(args) > 3){ 
        tagline= args[4] 
    } else {
        # Temporary fix:
        tagline = sub(".*/", "", bedfile)
        tagline = paste("File:", sub("_seed.*", "", tagline))
    }
}

# Ensure directory exists:
imgdir = paste0(resultdir, '/img/')
cmd = paste("mkdir -p", resultdir, imgdir)
system(cmd)


namesfile = paste0(filepref, '_names.tsv')
epinames = scan(namesfile,'c')
# is.epi = length(grep("BSS",epinames)) == nrow(gwmat)
# is.epi = length(grep("BSS",epinames)) == 833

if (!is.epi){
    # Load modules (TODO: Load extracted full epigenomes?)
    commandArgs <- function(trailingOnly=TRUE){
        c(filepref, tagline, imgdir) }
    source(paste0(bindir, 'load_modules_data.R'))
}


# ---------------------
# Load enrichment info:
# ---------------------
clong = read.delim(paste0(resultdir, '/enr_stats_table.tsv'), header=T)
data.file = paste0(resultdir, '/hg_processed_enr.Rda')
load(data.file)



# Process/plot a single matrix:
fdr = '99%'
rmat = calist[[fdr]]
# rmat = cplist[[fdr]]
kuid = rownames(rmat)
rmat[is.na(rmat)] = 0
kept.uids = names(which(apply(rmat > 0, 1, sum) > 0))
gwmat = t(rmat[kept.uids,])
rownames(gwmat) = paste0('c', c(1:nrow(gwmat) - 1))
keep.cls = list()
keep.cls[[tagline]] = rownames(gwmat)

# Choose # gwas to show:
SHOWGWAS = length(kuid)
zmax = 20; zmin = 4
gwasmarg = sort(apply(gwmat, 2, max), decreasing=T)
keep.studies = names(head(gwasmarg, SHOWGWAS))
# Order the top studies:
r2 = reord(t(gwmat[, keep.studies]) > zmin, 'Jaccard')
studyord = rownames(r2)
r3 = reord(gwmat[, keep.studies] > zmin, 'eJaccard')
roword = rownames(r3)
finalmat = gwmat[roword, studyord]
gmat = finalmat

# ---------------------------------
# Plot the gwas alone, no metadata:
# ---------------------------------
plot.alone = FALSE
if (plot.alone){
    # Diagonalize, ignore all <2
    tmat = gwmat
    tmat[tmat < zmin] = 0
    ll = diag.mat(tmat[roword,keep.studies])
    tmp = ll[[1]]
    cto = ll[[3]] # For breaks
    tmp[tmp > zmax] <- zmax

    # TODO: adjust height and width to SHOWGWAS and nrow(gmat)
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
}




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




if (is.epi){
















} else {
    # PLOT AS CLUSTERS:


    # Set groupdist parameters
    rnorder = rngroup # ord.gwas
    acutnamed = acutgroup.nam
    set = tagline 
    lablist = lablist.group
    dist.order = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, subset=TRUE, calc.only=TRUE)
    gll = plot.gwas(finalmat, set, dist.order, gwasord=studyord, calc.only=TRUE)
    gwnam = gll[[1]]
    gcut = gll[[2]]


    fdr.tag = paste0(100 - as.numeric(sub("%","",fdr)), 'pct')
    # Make plot:
    basic.plot = paste0(imgpref, "gwas_flathg_top", SHOWGWAS, "_", fdr.tag, "_groupord.png")
    png(basic.plot,res=450,units='in',width=12,height=17)
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
}














print(paste("[STATUS] Finished plotting for", bedfile))
