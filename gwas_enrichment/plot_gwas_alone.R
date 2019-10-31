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
source(paste0(bindir, 'load_metadata.R'))
today <- format(Sys.time(), "%m%d%y")
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

# Defaults:
gwasfile = 'observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_5000_enrich.tsv'
# gwasfile = 'observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_1000_enrich.tsv'
# gwasfile = 'observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_0_enrich.tsv'
extension = 5000
filepref = 'cls_merge2_wH3K27ac100_raw'
tagline = 'ChromHMM Enhancers (on Epigenomes)'
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

# Prefix / plot directory:
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'gwas_',filepref,'_e', extension, '_')
print(paste("[STATUS] Plotting under prefix", imgpref))

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
    # TODO: TRY REDUCE GWMAT IF THE SAME NAME?
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

    
    # ------------------------------------------------
    # Measure how many enrichments vs. old epigenomes:
    # ------------------------------------------------
    cutoff = 3
    passdf = filter(gwdf, logpval >= cutoff)  # Keep only passing cutoff - NOTE: try more stringent
    # Mapping:
    idmap = data.frame(id=epinames, cls=names(epinames))
    idmap = merge(idmap, meta[,c('id','GROUP','infoline','Project')])
    passdf = merge(passdf, idmap)

    oldepi = c('ENCODE 2012','Roadmap 2015')
    oldpassdf = filter(passdf, Project %in% oldepi)
    # Raw # of enrichments + %
    print(paste(nrow(oldpassdf),"to",nrow(passdf)))
    # New gwas studies (PMIDxGWAS) annotated # + %
    oldst = unique(oldpassdf$pmt)
    newst = unique(passdf$pmt)
    print(paste(length(oldst), "to", length(newst)))

    # TODO: All st:
    topenr = c()
    topnum = c()
    for (st in newst){
        stdf = filter(passdf, pmt == st)
        mrow = stdf[which.max(stdf$logpval),]
        sid = as.character(mrow$id)
        topenr = c(topenr, sid)
        topnum = c(topnum, mrow$logpval)
    }
    topdf = data.frame(id=topenr, num=topnum, pmt=newst)
    topdf = merge(topdf, idmap)

    tdf = aggregate(pmt ~ Project, topdf, length)
    tdf = tdf[order(tdf$pmt, decreasing=T),]
    tcols = colvals[['project']][tdf$Project]
    vals = tdf$pmt
    names(vals) = tdf$Project

    png(paste0(imgpref,'numtop_byproj_cutoff',cutoff,'.png'),res=450,units='in',width=6,height=2.5)
    par(mar=c(4.5,.5,1,3))
    par(yaxs='i')
    par(xaxs='i')
    barplot(tdf$pmt, width=.5, col=tcols, # names.arg=tdf$Project, 
            xlab='Number of GWAS', space=0.1, xlim = c(0, (round(max(tdf$pmt) / 50) + 1) * 50),
            main='# Strongest enrichment for each project', border=NA, horiz=T)
    abline(v=seq(0,max(tdf$pmt),50), col='white', lwd=2)
    txtoffset = diff(par()$usr[3:4]) / (2 * nrow(tdf))
    txtloc = seq(par()$usr[3] + txtoffset, par()$usr[4] - txtoffset, length.out=nrow(tdf))
    text(y=txtloc, x=tdf$pmt + 3, tdf$Project, 
         xpd=TRUE, srt=0, adj=0, cex=1)
    dev.off()

    # New studies:
    diffst = newst[!(newst %in% oldst)]
    # By strength of assoc:
    diffdf = filter(passdf, pmt %in% diffst)
    diffdf = diffdf[order(diffdf$logpval, decreasing=T),]
    head(diffdf[,c('trait','logpval','GROUP','infoline','Project')], 20)

    # Rarefaction-style curves (both types)
    rdf = passdf
    satord = c()
    satnum = c()
    while (nrow(rdf) != 0){
        ag = aggregate(pmt ~ id, rdf, length)
        mrow = ag[which.max(ag$pmt),]
        sid = as.character(mrow[[1]])
        npm = mrow[[2]]
        satord = c(satord, sid)
        satnum = c(satnum, npm)
        # Remove:
        pid = filter(rdf, id == sid)$pmt
        rdf = filter(rdf, !(pmt %in% pid))
    }
    satdf = data.frame(id=satord, num=satnum, cs=cumsum(satnum))
    satdf = merge(satdf, idmap)
    satdf = satdf[order(satdf$cs),]

    satcols = colvals[['project']][satdf$Project]
    satgroupcols = colvals[['group']][satdf$GROUP]
    png(paste0(imgpref,'saturation_all_cutoff',cutoff,'.png'),res=450,units='in',width=7,height=12)
    par(mar=c(4.5,4.5,1,1))
    plot(satdf$cs, rev(1:nrow(satdf)), pch=19, col=satcols, xlim=c(0, max(satdf$cs)), 
         ylab='Number of epigenomes', xlab='Cumulative number of GWAS')
    abline(v=seq(0,max(satdf$cs),50), col='darkgrey',lty='dotted')
    text(satdf$cs - 0.02 * diff(par()$usr[1:2]), rev(1:nrow(satdf)), labels=satdf$infoline, 
         srt=0, adj=1, cex=.5, col=satgroupcols)
    legend('bottom', names(colvals[['project']]), col=colvals[['project']], 
           pch=19,inset=.01, box.col=NA, title=expression(bold('Project')))
    legend('bottomleft', names(colvals[['group']]), col=colvals[['group']], 
           pch=19,inset=.01, box.col=NA, title=expression(paste(bold("Group"))))
    dev.off()

    # Rarefaction by project:
    projord = c("ENCODE 2012", "Roadmap 2015", 
                "Roadmap (Novel)", "ENCODE", "GGR")
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
        ag = aggregate(pmt ~ id, subdf, length)
        mrow = ag[which.max(ag$pmt),]
        sid = as.character(mrow[[1]])
        npm = mrow[[2]]
        satord = c(satord, sid)
        satnum = c(satnum, npm)
        # Remove:
        pid = filter(rdf, id == sid)$pmt
        rdf = filter(rdf, !(pmt %in% pid))
    }
    satdf = data.frame(id=satord, num=satnum, cs=cumsum(satnum))
    satdf = merge(satdf, idmap)
    satdf = satdf[order(satdf$cs),]

    # Which NEW TISSUES are the most informative?
    satcols = colvals[['project']][satdf$Project]
    satgroupcols = colvals[['group']][satdf$GROUP]
    png(paste0(imgpref,'saturation_byproj_cutoff',cutoff,'.png'),res=300,units='in',width=7,height=12)
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



}




