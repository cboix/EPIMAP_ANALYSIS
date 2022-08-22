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
# Or directory?
# filepref = 'cls_merge2_wH3K27ac100_raw'
# tagline = 'ChromHMM Enhancers (on Epigenomes)'
filepref = 'cls_merge2_wH3K27ac100_300'
tagline = 'ChromHMM Enhancers (on Modules)'
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

# Prefix / plot directory:
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'fimo_',filepref,'_', fset, '_')
print(paste("[STATUS] Plotting under prefix", imgpref))

# --------------------------------
# Add gwas enrichment information:
# --------------------------------
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
motdf$rmotif = redmot
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

# Reorder + plot motifs alone:
method='euclidean'
r2 = reord(t(holdmat),method)
r3 = reord(holdmat, method)
epiord = rownames(r2)
motord = rownames(r3)
finalmat = t(smat[motord, epiord])

# Diagonalize, ignore all <2
tmat = abs(smat[,epiord])
ll = diag.mat(t(tmat))
tmp = ll[[1]]
cto = ll[[3]] # For breaks
diagord = colnames(tmp)

# Plot gwas alone 
png(paste0(imgpref,'alone.png'),res=450,units='in',width=10,height=15)
par(mar=c(0.25, 6, .25, .25))
image(finalmat, axes=FALSE,col=rev(colrb), zlim=c(-tcut, tcut))
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
rows = (seq(0, ncol(finalmat), 10) - .5) / (ncol(finalmat) - 1)
cols = (seq(0, nrow(finalmat), 20) - .5) / (nrow(finalmat) - 1)
abline(h=rows,lty='dotted',lw=1, col='darkgrey')
abline(v=cols,lty='dotted',lw=1, col='darkgrey')
text(y=seq(0,1, length.out=length(motord)),
     x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=motord, srt=0, adj=1, xpd=TRUE,cex=.25)
dev.off()


png(paste0(imgpref,'alone_diag.png'),res=450,units='in',width=10,height=17)
par(mar=c(0.25, 6, .25, .25))
image(finalmat[,diagord], axes=FALSE,col=rev(colrb), zlim=c(-tcut, tcut))
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
# abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
rows = (seq(0, ncol(finalmat), 10) - .5) / (ncol(finalmat) - 1)
cols = (seq(0, nrow(finalmat), 20) - .5) / (nrow(finalmat) - 1)
abline(h=rows,lty='dotted',lw=1, col='darkgrey')
abline(v=cols,lty='dotted',lw=1, col='darkgrey')
text(y=seq(0,1, length.out=ncol(finalmat)),
     x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=diagord, srt=0, adj=1, xpd=TRUE,cex=.25)
dev.off()


# ------------------------------
# After plot alone, try matched:
# ------------------------------
if (is.epi){
    epimat = finalmat
    rawmat = t(holdmat)

    # Order
    method = 'euclidean'
    dt = dist(rawmat, method)
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

    # Set parameters
    rnorder = rnmat 
    acutnamed = acutmat.nam
    set = tagline 
    lablist = lablist.mat

    png(paste0(imgpref,'epi.png'),res=450,units='in',width=15,height=17)
    hbreaks = calc.breaks.acut(acutnamed)
    layout(matrix(c(1,2),1,2), widths=c(1.75,10), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    par(mar=c(5, 6, 2, 0.25))
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
    par(mar=c(5, 0.25, 2, 0.25))
    image(t(finalmat[rnorder, motord]), axes=FALSE,col=rev(colrb), zlim=c(-tcut, tcut))
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    rows = (seq(0, nrow(finalmat), 20) - .5) / (nrow(finalmat) - 1)
    abline(h=hbreaks,lty='dotted',lw=1, col='darkgrey')
    abline(v=rows,lty='dotted',lw=1, col='darkgrey')
    text(x=seq(0,1, length.out=length(motord)),
         y=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), 
         labels=motord, srt=90, adj=1, xpd=TRUE,cex=.25)
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

    png(paste0(imgpref,'epi_groupord.png'),res=450,units='in',width=15,height=17)
    hbreaks = calc.breaks.acut(acutnamed)
    layout(matrix(c(1,2),1,2), widths=c(1.75,10), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    par(mar=c(5, 6, 2, 0.25))
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
    par(mar=c(5, 0.25, 2, 0.25))
    image(t(finalmat[rnorder, motord]), axes=FALSE,col=rev(colrb), zlim=c(-tcut, tcut))
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    rows = (seq(0, nrow(finalmat), 20) - .5) / (nrow(finalmat) - 1)
    abline(h=hbreaks,lty='dotted',lw=1, col='darkgrey')
    abline(v=rows,lty='dotted',lw=1, col='darkgrey')
    text(x=seq(0,1, length.out=length(motord)),
         y=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), 
         labels=motord, srt=90, adj=1, xpd=TRUE,cex=.25)
    mtext(tagline, side=3, line=0.25)
    dev.off()


    # Diagonalize, ignore all <2
    tmat = abs(rawmat[rnorder,])
    ll = diag.mat(tmat)
    tmp = ll[[1]]
    cto = ll[[3]] # For breaks
    diagord = colnames(tmp)
    vcut =  c(cto[cto ==0], acutnamed[cto])
    vbreaks = calc.breaks.acut(vcut)

    png(paste0(imgpref,'epi_groupord_diag.png'),res=450,units='in',width=15,height=17)
    hbreaks = calc.breaks.acut(acutnamed)
    layout(matrix(c(1,2),1,2), widths=c(1.75,10), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    par(mar=c(5, 6, 2, 0.25))
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
    par(mar=c(5, 0.25, 2, 0.25))
    image(t(finalmat[rnorder, diagord]), axes=FALSE,col=rev(colrb), zlim=c(-tcut, tcut))
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    rows = (seq(0, nrow(finalmat), 20) - .5) / (nrow(finalmat) - 1)
    abline(h=hbreaks,lty='dotted',lw=1, col='darkgrey')
    abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
    text(x=seq(0,1, length.out=length(diagord)),
         y=par()$usr[3]-0.005*(par()$usr[4]-par()$usr[3]), 
         labels=diagord, srt=90, adj=1, xpd=TRUE,cex=.25)
    mtext(tagline, side=3, line=0.25)
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

    # ------------------------------------
    # Flip and add rectangles on diagonal:
    png(paste0(imgpref,'epi_groupord_diag_flip.png'),res=450,units='in',width=16,height=17)
    layout(matrix(c(1,2),2,1), heights=c(1.1,10), TRUE)
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
    par(mar=c(.25, 15, .25, 0.25))
    image(finalmat[rnorder, rev(diagord)], axes=FALSE,col=rev(colrb), zlim=c(-tcut, tcut))
    # image(epimat[rnorder, rev(diagord)], axes=FALSE,col=colred, zlim=c(zmin, zmax))
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    hbreaks = calc.breaks.acut(acutnamed)
    vbreaks = calc.breaks.acut(vcut)
    abline(v=hbreaks,lty='dotted',lw=1, col='darkgrey')
    abline(h=par()$usr[2] - vbreaks,lty='dotted',lw=1, col='darkgrey')
    text(y=seq(0,1, length.out=length(diagord)),
         x=par()$usr[1]-0.005*(par()$usr[2]-par()$usr[1]), 
         labels=rev(diagord), srt=0, adj=1, xpd=TRUE,cex=.25)
    # Add rectangles
    rll = calc.breaks.rect(hcls=rev(vcut), vcls=acutnamed, colset)
    rectdf = rll[[1]]
    vccols = rll[[2]]
    # Add grey - ubq rectangle:
    rectdf = rbind(c(x1=par()$usr[1], x2=par()$usr[2],
                     y1=rectdf$y2[1], y2=par()$usr[4]), rectdf)
    vccols = c('grey',vccols)
    rect(xleft=rectdf$x1, xright=rectdf$x2,
         ybottom=rectdf$y1, ytop=rectdf$y2,
         border=vccols, lwd=1)
    par(xpd=NA)
    rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
         xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
         border='white', col=vccols, lwd=.25)
    # rect(xleft=rectdf$x1, xright=rectdf$x2, ybottom=par()$usr[3], 
    #      ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
    #      border='white', col=vccols, lwd=.25)
    dev.off()

}

