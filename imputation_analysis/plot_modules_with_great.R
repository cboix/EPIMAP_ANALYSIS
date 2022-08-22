#!/usr/bin/R
# -----------------------------------------------------
# Plot modules from CHMM/other clustering
# With GREAT
# TODO: 
# 1. Plot the centers matrix with great
# 2. Plot stats on great - amt of enrichment/depletion
# 3. Associate great with cell types by enrichment
# 4. Plot and order/cluster the adjacency matrix
#       for the epigenomes by great co-enrichment
# -----------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
# source(paste0(bindir, 'great_enrichment/general_WM_functions.R'))
source(paste0(bindir, 'great_enrichment/auxfunctions_GREAT.R'))

# Defaults:
resultsdir = './'
filepref = 'cls_merge2_wH3K27ac100_300'
tagline = 'ChromHMM Enhancers'
# resultsdir = './prom_cls/'
# filepref = 'prom_cls/cls_merge2_wH3K27ac100_300'
# tagline = 'ChromHMM Promoters'
imgdir = paste0(img, "clusters/") 

# Arguments:
args=(commandArgs(TRUE))
if (length(args) == 0) {
    stop("No arguments supplied: Need clusters filename")
} else {        
    resultsdir = args[1]
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
subimgdir = paste0(imgdir, 'go_figures/')
cmd = paste('mkdir -p', imgdir, subimgdir)
system(cmd)
imgpref = paste0(imgdir, 'clusters_great_', sub("/","_", filepref), '_')

# ------------------------------
# Load ontology enrichment data:
# ------------------------------
parseddir <- paste(resultsdir, "parsed_results", sep="/");
load(paste(parseddir, "/res_all_per_onto.RData", sep=""))
ontos <- names(res_all)  # Ontology names
# Keep only some ontologies (core or extended):
coreonto = c('GO Biological Process', 'GO Molecular Function', 
             'GO Cellular Component', 'Human Phenotype')
goonto= c('GO Biological Process', 'GO Molecular Function', 'GO Cellular Component')
extonto = c(coreonto, 'Mouse Phenotype', 'Disease Ontology')
if (length(goonto) == 1){
    res_kept = res_all[[goonto]]
} else {
    res_kept = sapply(goonto, function(x){res_all[[x]]})
}
# res_kept = sapply(coreonto, function(x){res_all[[x]]})
# res_kept = sapply(extonto, function(x){res_all[[x]]})

# Abbreviate the names of the rows:
abbrev <- sapply(names(res_kept), function(x) {
                     paste(sapply(strsplit(x, ""), function(y) { 
                                      sel <- which(y %in% LETTERS);
                                      sel <- sort(unique(c(sel, sel+1)));
                                      sel <- sel[y[sel] %in% LETTERS | y[sel] %in% letters]
                                      y[sel] }), collapse="") })
ontonames = sapply(names(res_kept), function(x){ paste(abbrev[x], rownames(res_kept[[x]]), sep=" - ")})
rownames_filt <- as.character(unlist(ontonames))

# ----------------------
# Create ontology matrix
# ----------------------
if (length(goonto) == 1){
    resmat <- as.matrix(res_kept)
} else {
    res_kept_concat <- do.call(rbind, res_kept)
    resmat <- as.matrix(res_kept_concat)
}
colnames(resmat) = paste0('c', colnames(resmat))
# Filter resmat:
cutoff = 2.5
resmarg = apply(resmat,1,max)
keepres = which(resmarg > cutoff)
resmat = resmat[keepres,]

write.table(resmat,paste0(parseddir, '/', filepref, '_kept_go.tsv'), col.names=T, row.names=T, sep="\t", quote=F)

# Change names if they correspond to epigenomes:
namesfile = paste0(filepref, '_names.tsv')
epinames = scan(namesfile,'c')
is.epi = length(grep("BSS",epinames)) == ncol(resmat)
if (is.epi){
    print("Is enrichment on epigenomes, changing rownames epigenomes")
    clsn = paste0('c',1:length(epinames) - 1)
    names(epinames) = clsn
    colnames(resmat) = epinames[colnames(resmat)]
}

# -------------------------------------------
# Order the complete matrix by binary/jaccard
# distance of significant results.
# -------------------------------------------
noreorder = FALSE
method = 'ejaccard'
resmat_row_ord <- get_optimal_ordering(dist(resmat > 2, method=method))
if (noreorder){
    resmat_col_ord <- 1:ncol(resmat)
} else {
    resmat_col_ord <- get_optimal_ordering(dist(t(resmat > 2), method=method));
}
resmat <- resmat[resmat_row_ord,]

# Save reordered + selected rownames and get wordlets:
rownames_filt <- rownames_filt[resmat_row_ord]
save(resmat, rownames_filt, file=paste0(resultsdir, "resmat_rownames_filt.RData"))


# --------------------------------
# Set up plot with group ordering:
# --------------------------------
# Variables for plot order:
rnorder = rngroup
acutnamed = acutgroup.nam
set = tagline 
lablist = lablist.group

# Get cluster order and diagonalize GO terms:
dist.order = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, calc.only=TRUE)
clsord = dist.order[[1]]
clscuts = dist.order[[3]]

# Diagonalize, ignore all <zmin
tmin = cutoff
tmat = resmat[,clsord]
tmat[tmat < tmin] = 0
ll = diag.mat(t(tmat))
tmp = ll[[1]]
cto = ll[[3]] # For breaks
diagord = colnames(tmp)
dresmat = resmat[diagord,clsord]
# Use cto to get breaks:
vcut2 =  c(cto[cto ==0], clscuts[cto])
hbreaks = calc.breaks.acut(vcut2)

# Get appropriate wordlets for diagonal matrix:
NBIN = 10
wordlets <- get_summary_terms(diagord, NBIN)
wordlets_pos <- seq(1, length(diagord), length.out=length(wordlets)+1);
wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;

plot.great <-  function(gmat, set, ordll, hbreaks=NULL, colp=TRUE,
                        zmin=0.5, zmax=5, labeled=TRUE, ablwd=1){
    ord = ordll[[1]]
    vbreaks = ordll[[2]]
    subord = ord[ord %in% colnames(gmat)]
    gmat = t(gmat[,subord])
    # Threshold for plotting:
    gmat[gmat > zmax] <- zmax
    gmat[gmat < zmin] <- 0
    if (colp){
        pmat = gmat
        pmat[pmat < 2] = 1
        pmat[pmat < 3 & pmat >= 2] = 2
        pmat[pmat < 4 & pmat >= 3] = 3
        pmat[pmat >= 4] = 4
        rcols = c('white','yellow','orange','red')
        image(pmat, axes=FALSE,col=rcols, zlim=c(zmin, zmax), useRaster=TRUE)
    } else {
        image(gmat, axes=FALSE,col=col, zlim=c(zmin, zmax), useRaster=TRUE)
    }
    # Lines:
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    rows = (seq(0, ncol(gmat), 200) - .5) / (ncol(gmat) - 1)
    abline(v=vbreaks,lty='dotted',lw=ablwd, col='darkgrey')
    if (!is.null(hbreaks)){
        abline(h=hbreaks,lty='dotted',lw=ablwd, col='darkgrey')
    }
}


png(paste0(imgpref,'n',NBIN,'_groupdist.png'),res=450,units='in',width=12,height=17)
ratio = 1
layout(matrix(c(1:12),4,3), heights=c(.6,8,8 * ratio, 1), widths=c(1.25,8,.5), TRUE)
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
image(t(dresmat[,clsord]), col='white', axes=F, useRaster=TRUE)
text(y=wordlets_pos / nrow(dresmat),
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=wordlets, srt=0, adj=1, xpd=TRUE,cex=.3)
# Add labels
plot.new()
par(mar=c(0.25, 0.25, 2, 0.25))
meta.image(enrichmat[clsord,5:1], colvals=colvals, cex=0, horiz=F, useRaster=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
vbreaks = dist.order[[set]][[2]]
abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
mtext(set, side=3, cex=1.3)
# Plot clusters and counts
set = tagline
par(mar=c(0.25, 0.25, 0, 0.25))
dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, title=FALSE)
abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
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
par(mar=c(0.25, 0.25, 0, 0.25))
plot.great(dresmat, set, dist.order[[set]], hbreaks=hbreaks)
# Add rectangles:
rll = calc.breaks.rect(hcls=vcut2, vcls = dist.order[[set]][[3]], colset)
rectdf = rll[[1]]
vccols = rll[[2]]
rect(xleft=rectdf$x1, xright=rectdf$x2,
     ybottom=rectdf$y1, ytop=rectdf$y2,
     border=vccols, lwd=1)
par(xpd=NA)
rect(xleft=rectdf$x1, xright=rectdf$x2, ybottom=par()$usr[3], 
     ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=vccols, lwd=.25)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
par(xpd=FALSE)
# Add counts:
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



# Keeping only "specific GO terms"
tmin = cutoff
ind = which(apply(resmat >= 2, 1, sum) < 30 & apply(resmat, 1, max) >= 4)
tmat = resmat[ind, clsord]
tmat[tmat < tmin] = 0

ll = diag.mat2(t(tmat), ratio=.1)
tmp = ll[[1]]
cto = ll[[3]] # For breaks
diagord = colnames(tmp)
dresmat = resmat[diagord,clsord]
dr2resmat = resmat[diagord,clsord]
vcut2 =  c(cto[cto ==0], clscuts[cto])
hbreaks = calc.breaks.acut(vcut2)

# Get appropriate wordlets for diagonal matrix:
NBIN = 10
diagord2 = diagord
diagord2[sapply(diagord2, nchar) > 60] = ''
wordlets <- get_summary_terms(diagord2, NBIN)
wordlets_pos <- seq(1, length(diagord), length.out=length(wordlets)+1);
wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;

# Color wordlets by category:
wind = sapply(wordlets, function(x){which(diagord == x)})
ctw = vcut2[wind]
wcol = c('grey', colset)[ctw + 1]

png(paste0(imgpref,'n',NBIN,'_groupdist_specific_onlygo.png'),res=450,units='in',width=7.75,height=3.5)
# pdf(paste0(imgpref,'n',NBIN,'_groupdist_specific_onlygo.pdf'), width=7.75, height=3.5)
layout(matrix(c(1:2),1,2), widths=c(2,8), TRUE)
par(xpd=FALSE)
par(yaxs="i")
par(xaxs="i")
sp=0.15
dist.order = list()
# GWAS labels:
par(mar=c(sp, 5, sp, sp))
image(t(dresmat[,clsord]), col='white', axes=F, useRaster=TRUE)
# text(y=wordlets_pos / nrow(dresmat),
text(y=seq(0+1e-3,1-1e-3, length.out=length(wordlets)),
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=wordlets, srt=0, adj=1, xpd=TRUE,cex=.3, col=wcol)
# Load matrix:
dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, title=FALSE, calc.only=TRUE)
# Plot go enrichments:
par(mar=rep(sp,4))
plot.great(dresmat, set, dist.order[[set]], hbreaks=hbreaks, ablwd=.5)
rll = calc.breaks.rect(hcls=vcut2, vcls = dist.order[[set]][[3]], colset)
rectdf = rll[[1]]
vccols = rll[[2]]
rect(xleft=rectdf$x1, xright=rectdf$x2,
     ybottom=rectdf$y1, ytop=rectdf$y2,
     border=vccols, lwd=.75)
par(xpd=NA)
# Cls horizontal boxes:
rll = calc.breaks.rect(hcls=acutnamed, vcls=dist.order[[set]][[3]], colset)
clsrectdf = rll[[1]]
clsvccols = rll[[2]]
# Add grey - ubq rectangle:
clsrectdf = rbind(c(x1=par()$usr[1], x2=clsrectdf$x1[1],
                 y1=par()$usr[3], y2=par()$usr[4]), clsrectdf)
clsvccols = c('grey',clsvccols)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[4], 
     ytop=par()$usr[4]+0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[3], 
     ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
par(xpd=FALSE)
# TODO: Add legend. (or on motifs figure)
dev.off()


# ------------------
# Plot all go terms:
# ------------------
wordlets = diagord
wordlets_pos <- seq(0, 1, length.out=length(diagord))

# Color wordlets by category:
wind = sapply(wordlets, function(x){which(diagord == x)})
ctw = vcut2[wind]
wcol = c('grey', colset)[ctw + 1]

midpt = round(nrow(dresmat)/2)
botid = 1:midpt
topid = (midpt + 1):nrow(dresmat)

plot_all_go_splitid = function(pltid, sp=0.1, txtcex=.15){
    layout(matrix(c(1:2),1,2), widths=c(2.5,8), TRUE)
    par(xpd=FALSE)
    par(yaxs="i")
    par(xaxs="i")
    dist.order = list()
    # GWAS labels:
    par(mar=c(sp, 5, sp, sp))
    image(t(dresmat[pltid,clsord]), col='white', axes=F, useRaster=TRUE)
    # text(y=wordlets_pos / nrow(dresmat),
    text(y=seq(0,1, length.out=length(diagord[pltid])),
         x=par()$usr[2]-0.0001*(par()$usr[2]-par()$usr[1]), 
         labels=diagord[pltid], srt=0, adj=1, xpd=TRUE,cex=txtcex, col=wcol[pltid])
    # Load matrix:
    dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, title=FALSE, calc.only=TRUE)
    # Plot go enrichments:
    par(mar=rep(sp,4))
    hbreaks.local = calc.breaks.acut(vcut2[pltid])
    plot.great(dresmat[pltid,], set, dist.order[[set]], hbreaks=hbreaks.local, ablwd=.5)
    rll = calc.breaks.rect(hcls=vcut2[pltid], vcls = dist.order[[set]][[3]], colset)
    rectdf = rll[[1]]
    vccols = rll[[2]]
    rect(xleft=rectdf$x1, xright=rectdf$x2,
         ybottom=rectdf$y1, ytop=rectdf$y2,
         border=vccols, lwd=.75)
    par(xpd=NA)
    # Cls horizontal boxes:
    rll = calc.breaks.rect(hcls=acutnamed, vcls=dist.order[[set]][[3]], colset)
    clsrectdf = rll[[1]]
    clsvccols = rll[[2]]
    # Add grey - ubq rectangle:
    clsrectdf = rbind(c(x1=par()$usr[1], x2=clsrectdf$x1[1],
                        y1=par()$usr[3], y2=par()$usr[4]), clsrectdf)
    clsvccols = c('grey',clsvccols)
    rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[4], 
         ytop=par()$usr[4]+0.004*(par()$usr[4]-par()$usr[3]), 
         border='white', col=clsvccols, lwd=.25)
    rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[3], 
         ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
         border='white', col=clsvccols, lwd=.25)
    rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
         xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
         border='white', col=vccols, lwd=.25)
    par(xpd=FALSE)
    # TODO: Add legend. (or on motifs figure)
}

pdf(paste0(imgpref,'all_groupdist_specific_onlygo_tophalf.pdf'), width=8.25, height=10)
plot_all_go_splitid(topid)
dev.off()

pdf(paste0(imgpref,'all_groupdist_specific_onlygo_bothalf.pdf'), width=8.25, height=10)
plot_all_go_splitid(botid)
dev.off()


# png(paste0(imgpref,'all_groupdist_specific_onlygo.png'),res=450,units='in',width=7.75,height=3.5)
pdf(paste0(imgpref,'all_groupdist_specific_onlygo.pdf'), width=8.25, height=25)
dev.off()





# Plot all, for supplement:
# Keeping only "specific GO terms"
tmin = cutoff
ind = which(apply(resmat >= 2, 1, sum) < 30 & apply(resmat, 1, max) >= 4)
tmat = resmat[ind, clsord]
tmat[tmat < tmin] = 0

ll = diag.mat2(t(tmat), ratio=.1)
tmp = ll[[1]]
cto = ll[[3]] # For breaks
diagord = colnames(tmp)
dresmat = resmat[diagord,clsord]
dr2resmat = resmat[diagord,clsord]
vcut2 =  c(cto[cto ==0], clscuts[cto])
hbreaks = calc.breaks.acut(vcut2)

# Get appropriate wordlets for diagonal matrix:
NBIN = 1
# diagord2 = diagord
# diagord2[sapply(diagord2, nchar) > 40] = ''
# wordlets <- get_summary_terms(diagord2, NBIN)
wordlets = diagord
wordlets_pos <- seq(1, length(diagord), length.out=length(wordlets)+1);
wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;

# Color wordlets by category:
wind = sapply(wordlets, function(x){which(diagord == x)})
ctw = vcut2[wind]
wcol = c('grey', colset)[ctw + 1]


png(paste0(imgpref,'n',NBIN,'_groupdist_specific_onlygo.png'),res=450,units='in',width=7.75,height=30)
layout(matrix(c(1:2),1,2), widths=c(4,8), TRUE)
par(xpd=FALSE)
par(yaxs="i")
par(xaxs="i")
sp=0.15
dist.order = list()
# GWAS labels:
par(mar=c(sp, 5, sp, sp))
image(t(dresmat[,clsord]), col='white', axes=F, useRaster=TRUE)
text(y=seq(0,1, length.out=length(wordlets)),
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=wordlets, srt=0, adj=1, xpd=TRUE,cex=.2, col=wcol)
# Load matrix:
dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, title=FALSE, calc.only=TRUE)
# Plot go enrichments:
par(mar=rep(sp,4))
plot.great(dresmat, set, dist.order[[set]], hbreaks=hbreaks, ablwd=.5)
rll = calc.breaks.rect(hcls=vcut2, vcls = dist.order[[set]][[3]], colset)
rectdf = rll[[1]]
vccols = rll[[2]]
rect(xleft=rectdf$x1, xright=rectdf$x2,
     ybottom=rectdf$y1, ytop=rectdf$y2,
     border=vccols, lwd=.75)
par(xpd=NA)
# Cls horizontal boxes:
rll = calc.breaks.rect(hcls=acutnamed, vcls=dist.order[[set]][[3]], colset)
clsrectdf = rll[[1]]
clsvccols = rll[[2]]
# Add grey - ubq rectangle:
clsrectdf = rbind(c(x1=par()$usr[1], x2=clsrectdf$x1[1],
                 y1=par()$usr[3], y2=par()$usr[4]), clsrectdf)
clsvccols = c('grey',clsvccols)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[4], 
     ytop=par()$usr[4]+0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[3], 
     ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
par(xpd=FALSE)
# TODO: Add legend. (or on motifs figure)
dev.off()


# ----------------------------------------
# Same small figure with better labelling:
# ----------------------------------------
tmin = cutoff
ind = which(apply(resmat >= 2, 1, sum) < 30 & apply(resmat, 1, max) >= 4)
tmat = resmat[ind, clsord]
tmat[tmat < tmin] = 0

ll = diag.mat2(t(tmat), ratio=.1)
tmp = ll[[1]]
cto = ll[[3]] # For breaks
diagord = colnames(tmp)
dresmat = resmat[diagord,clsord]
dr2resmat = resmat[diagord,clsord]
vcut2 =  c(cto[cto ==0], clscuts[cto])
hbreaks = calc.breaks.acut(vcut2)

# Count how many:
NTERMS = 80 # total number terms
termcex = 0.3
# NTERMS = 50 # total number terms
# termcex = 0.4
NREP = NTERMS - length(unique(vcut2)) # 1 per at least.
rdf = aggregate(count ~ cls, data.frame(cls = factor(vcut2, levels =unique(vcut2)), count=1), sum)
rdf$rep = round(rdf$count / (length(vcut2) / NREP)) + 1
rdf$cls = as.numeric(as.character(rdf$cls))
term_words <- strsplit(diagord, "[ _,.]");
tab_all <- sort(table(unlist(term_words)));
terms = c()
wcls = c()
for (i in 1:nrow(rdf)){
    wind = which(vcut2 == rdf$cls[i])
    subrn = diagord[wind]
    # if (length(subrn) > 5){ filt=TRUE } else { filt=FALSE }
    if (length(subrn) > 3){ filt=TRUE } else { filt=FALSE }
    sterms = get_summary_terms(subrn, binsize=round(length(subrn) / rdf$rep[i]), 
                               filtering=filt, tab_all=tab_all)
    terms = c(terms, sterms)
    wcls = c(wcls, rep(rdf$cls[i], length(sterms)))
}

# Color wordlets by category:
wordlets = terms
wcol = c('grey', colset)[wcls + 1]
wind = sapply(wordlets, function(x){which(diagord == x)})


# png(paste0(imgpref,'n',NTERMS,'_groupdist_specific_onlygo_repr.png'),res=450,units='in',width=7.75,height=3.5)
pdf(paste0(imgpref,'n',NTERMS,'_groupdist_specific_onlygo_repr.pdf'), width=8,height=3.5)
if (termcex > .3){
    layout(matrix(c(1:2),1,2), widths=c(2.5,8), TRUE)
} else {
    layout(matrix(c(1:2),1,2), widths=c(2.25,8), TRUE)
}
par(xpd=FALSE)
par(yaxs="i")
par(xaxs="i")
sp=0.15
dist.order = list()
par(mar=c(sp, 6, sp, 0))
# Add the GO labels:
image(t(dresmat[,clsord]), col='white', axes=F, useRaster=TRUE)
par(xpd=TRUE)
rx=c(0.005,0.05, 0.2,0.225)
xl = wordlets
xc = wcol
xx = seq(0+1e-3, 1-1e-3, length.out=length(xl))
xw = wind / (length(diagord) +1)
x = par()$usr[2]-rx*(diff(par()$usr[1:2]))
text(y=xx, x=x[4], labels=xl,
     srt=0, adj=1, xpd=TRUE, cex=termcex, col=wcol)
segments(x0=x[3], y0=xx, x1=x[2], y1=xw, col=wcol, lwd=.5)
segments(x0=x[1], y0=xw, x1=x[2], y1=xw, col=wcol, lwd=.5)
par(xpd=FALSE)
# Load matrix:
dist.order[[set]] = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, title=FALSE, calc.only=TRUE)
# Plot go enrichments:
par(mar=c(sp,sp,sp,sp))
plot.great(dresmat, set, dist.order[[set]], hbreaks=hbreaks, ablwd=.5)
rll = calc.breaks.rect(hcls=vcut2, vcls = dist.order[[set]][[3]], colset)
rectdf = rll[[1]]
vccols = rll[[2]]
rect(xleft=rectdf$x1, xright=rectdf$x2,
     ybottom=rectdf$y1, ytop=rectdf$y2,
     border=vccols, lwd=.75)
par(xpd=NA)
# Cls horizontal boxes:
rll = calc.breaks.rect(hcls=acutnamed, vcls=dist.order[[set]][[3]], colset)
clsrectdf = rll[[1]]
clsvccols = rll[[2]]
# Add grey - ubq rectangle:
clsrectdf = rbind(c(x1=par()$usr[1], x2=clsrectdf$x1[1],
                 y1=par()$usr[3], y2=par()$usr[4]), clsrectdf)
clsvccols = c('grey',clsvccols)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[4], 
     ytop=par()$usr[4]+0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
rect(xleft=clsrectdf$x1, xright=clsrectdf$x2, ybottom=par()$usr[3], 
     ytop=par()$usr[3]-0.004*(par()$usr[4]-par()$usr[3]), 
     border='white', col=clsvccols, lwd=.25)
rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
     xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
     border='white', col=vccols, lwd=.25)
par(xpd=FALSE)
dev.off()

