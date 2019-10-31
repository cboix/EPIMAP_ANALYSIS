#!/usr/bin/R
# -----------------------------------------------------
# Plot modules from CHMM/other clustering with GREAT
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
source(paste0(bindir, 'great_enrichment/auxfunctions_GREAT.R'))
today <- format(Sys.time(), "%m%d%y")
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

# Defaults:
resultsdir = './'
filepref = 'cls_merge2_wH3K27ac100_raw'
tagline = 'ChromHMM Enhancers (on Epigenomes)'
imgdir = paste0(img, "clusters/") 

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need clusters filename")
} else {        
    resultsdir = args[1]
    filepref = args[2]
    tagline = args[3]
    if (length(args) > 3){ 
        imgdir = args[4] 
    }
}

# Prefix:
subimgdir = paste0(imgdir, 'go_figures/')
cmd = paste('mkdir -p', imgdir, subimgdir)
system(cmd)
imgpref = paste0(imgdir, 'great_', filepref,'_')

# ------------------------------
# Load ontology enrichment data:
# ------------------------------
parseddir <- paste(resultsdir, "parsed_results", sep="/");
load(paste(parseddir, "/res_all_per_onto.RData", sep=""))
ontos <- names(res_all)  # Ontology names
# Keep only some ontologies (core or extended):
coreonto = c('GO Biological Process', 'GO Molecular Function', 
             'GO Cellular Component', 'Human Phenotype')
extonto = c(coreonto, 'Mouse Phenotype', 'Disease Ontology')
res_kept = sapply(coreonto, function(x){res_all[[x]]})
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
res_kept_concat <- do.call(rbind, res_kept)
resmat <- as.matrix(res_kept_concat)
colnames(resmat) = paste0('c', colnames(resmat))

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

# ---------------------------------------------
# Get wordlets:
# NOTE: Wordlets must change with diff ordering
# ---------------------------------------------
NBIN = 10
wordlets <- get_summary_terms(rownames(resmat), NBIN)
wordlets_pos <- seq(1, nrow(resmat), length.out=length(wordlets)+1);
wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;

# ------------------------------
# Plot great terms on their own:
# ------------------------------
# Threshold for plotting:
plot.great = function(pmat, zmax=5, zmin=0.5){
    pmat[pmat > zmax] <- zmax
    pmat[pmat < zmin] <- 0
    pmat[pmat < 2] = 1
    pmat[pmat < 3 & pmat >= 2] = 2
    pmat[pmat < 4 & pmat >= 3] = 3
    pmat[pmat >= 4] = 4
    rcols = c('white','yellow','orange','red')
    image(t(pmat), axes=FALSE,col=rcols, zlim=c(zmin, zmax))
    box(lwd=0.5)
}

rows = (seq(0, ncol(resmat), 50) - .5) / (ncol(resmat) - 1)
cols = (seq(0, nrow(resmat), 75) - .5) / (nrow(resmat) - 1)

png(paste0(imgpref,'n',NBIN,'_alone.png'),res=450,units='in',width=11,height=15)
layout(matrix(c(1),1,1), TRUE)
par(yaxs="i")
par(xaxs="i")
par(mar=c(0.25, 15, 2, 0.25))
plot.great(resmat[,resmat_col_ord])
# Lines:
abline(h=rows,lty='dotted',lw=1, col='darkgrey')
abline(v=cols,lty='dotted',lw=1, col='darkgrey')
text(y=wordlets_pos / nrow(resmat),
     x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=wordlets, srt=0, adj=1, xpd=TRUE,cex=.4)
mtext(tagline, side=3, line=0.5)
dev.off()


# Diagonalize, ignore all <2
tmin=2
tmat = resmat[,resmat_col_ord]
tmat[tmat < tmin] = 0
ll = diag.mat(t(tmat))
tmp = ll[[1]]
cto = ll[[3]] # For breaks
diagord = colnames(tmp)

# NOTE: Wordlets must change with diff ordering
dwordlets <- get_summary_terms(diagord, NBIN)
dwordlets_pos <- seq(1, nrow(resmat), length.out=length(dwordlets)+1);
dwordlets_pos <- dwordlets_pos[-1] - diff(dwordlets_pos)/2;

png(paste0(imgpref,'n',NBIN,'_alone_diag.png'),res=450,units='in',width=11,height=15)
layout(matrix(c(1),1,1), TRUE)
par(yaxs="i")
par(xaxs="i")
par(mar=c(0.25, 15, 2, 0.25))
plot.great(resmat[diagord,resmat_col_ord])
# Lines:
abline(h=rows,lty='dotted',lw=1, col='darkgrey')
abline(v=cols,lty='dotted',lw=1, col='darkgrey')
text(y=dwordlets_pos / nrow(resmat),
     x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=dwordlets, srt=0, adj=1, xpd=TRUE,cex=.4)
mtext(tagline, side=3, line=0.5)
dev.off()


# Plot for epigenomes:
if (is.epi){
    # Order
    method = 'ejaccard'
    dt = dist(t(resmat) > 2, method)
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
    # image(resmat[, rnmat], axes=FALSE,col=colred) 

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


    png(paste0(imgpref,'n',NBIN,'_epi.png'),res=450,units='in',width=15,height=17)
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
    plot.great(resmat[,rnmat])
    abline(h=rows,lty='dotted',lw=1, col='darkgrey')
    abline(v=cols,lty='dotted',lw=1, col='darkgrey')
    text(y=wordlets_pos / nrow(resmat),
         x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
         labels=wordlets, srt=0, adj=1, xpd=TRUE,cex=.5)
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

    # Diagonalize, ignore all <2
    tmat = resmat
    tmat[tmat < tmin] = 0
    ll = diag.mat(t(tmat[,rnorder]))
    ll = diag.mat(t(tmat[,rnorder]), ratio=0.25)
    tmp = ll[[1]]
    cto = ll[[3]] # For breaks
    diagord = rev(colnames(tmp))
    vcut =  c(cto[cto ==0], acutnamed[cto])
    vbreaks = calc.breaks.acut(rev(vcut))
    hbreaks = calc.breaks.acut(acutnamed)

    # NOTE: Wordlets must change with diff ordering
    NBIN = 8
    wordlets <- get_summary_terms(diagord, NBIN)
    wordlets_pos <- seq(1, nrow(resmat), length.out=length(wordlets)+1);
    wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;


    png(paste0(imgpref,'n',NBIN,'_epi_groupord_diag.png'),res=450,units='in',width=15,height=17)
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
    plot.great(resmat[diagord,rnorder])
    text(y=wordlets_pos / nrow(resmat),
         x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
         labels=wordlets, srt=0, adj=1, xpd=TRUE,cex=.5)
    abline(v=hbreaks,lty='dotted',lw=1, col='darkgrey')
    abline(h=vbreaks,lty='dotted',lw=1, col='darkgrey')
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
